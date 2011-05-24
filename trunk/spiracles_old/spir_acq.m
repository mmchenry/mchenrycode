function spir_acq

%% Set paths

%root = 'D:\PhD\research\spiracles\synchronisation';
root = 'C:\Users\Guest\Documents\spir';

% Syntax for filenames of video frame files
capFileName = 'Frame_$FI';

%% Find out experiment type

eType = questdlg('What kind of experiment do you want to analyze?',...
               'Experiment type','2 camera','camera + sensor','2 camera');

if isempty(eType)
    return
end

%% Instructions

disp(' ')
disp('Instructions __________________________________________________')
disp('   Control camera with Fire-i software');
disp('   Be sure Frame capture and trigger are clicked before recording');
disp('   Fire-i Settings:');
disp('   Frame capture - path should be set to root\video_tmp ')
disp('                 - Capture file format: TIF');
disp('                 - Number of frames to capture need to exceed')
disp('                    (frame rate) x (duration in sec)')

if strcmp(eType,'2 camera')
    disp(['                 - Capture file name: ' capFileName])
    disp(' ')
    disp('   cam0 -------------------')
end

disp('   Trigger - Polarity = High')
disp('           - Trigger Mode 15')
disp(' ')

%% Prompt for parameters 

% Prompt user for parameters
if strcmp(eType,'2 camera')
    prompt   = {'Individual number','Spiracle number','Recording number','Duration (s):',...
        'Frame rate (fps)','Body side of pike','Body side of marlin'};
    defaults = {'1','1','2','1','15','R','L'};
else
    prompt   = {'Individual number','Spiracle number','Recording number','Duration (s):',...
        'Frame rate (fps)','DAQ lower V','DAQ upper V'};
    defaults = {'1','1','2','1','15','-5','5'};
end

answer   = inputdlg(prompt,'Recording parameters',1,defaults);

% Grab individual and sequence numbers
indivNum = str2num(answer{1});
seqNum   = str2num(answer{3});

% Spiracle number
spirNum = str2num(answer{2});

% Recording duration (s)
duration = str2num(answer{4});

% Video recording rate (fps/Hz)
frame_rate = str2num(answer{5});

if strcmp(eType,'2 camera')
    % Side for each camera
    camSide0 = str2num(answer{6});
    camSide1 = str2num(answer{7});
    lowerV = nan;
    upperV = nan;
else
    % DAQ voltage range
    lowerV = str2num(answer{6});
    upperV = str2num(answer{7});
    camSide0 = ' ';
    camSide1 = ' ';
end

% Output sample rate (Hz) [triggers the camera]
sample_rate_out = 150;

% Input sample rate (Hz) [recording from optical sensor]
if strcmp(eType,'camera + sensor')
    sample_rate_in  = 10e2;
else
    sample_rate_in = [];
end

% Output voltage (V) [Note: do not exceed 5V!]
Vout = 3;

% Calculate number of frames
% if strcmp(eType,'2 camera')
%     numFrames = 2 * duration .* frame_rate;
% else
    numFrames = duration .* frame_rate;
%end

clear prompt defaults answer

%% Initialize input channel(s)

if strcmp(eType,'camera + sensor')
    AI     = analoginput('nidaq','Dev1');
    chanIn = addchannel(AI,0);
    set(AI,'SampleRate',sample_rate_in)
    
    ActualRate = get(AI,'SampleRate');
    set(AI,'SamplesPerTrigger',duration*ActualRate)
    set(AI,'TriggerType','Manual')
    blocksize = get(AI,'SamplesPerTrigger');
end

%% Initialize output channel

% Check output voltage
if Vout > 5
    beep;
    disp('Output voltage cannot exceed 5 V');
    delete(AI)
    clear AI
    return
end

AO = analogoutput('nidaq','Dev1');
chanOut = addchannel(AO,0);

set(AO,'SampleRate',sample_rate_out)
set(AO,'TriggerType','Manual')
ActualRate = get(AO,'SampleRate');

% Rezero the trigger channel & pause briefly to register
putsample(AO,0)
pause(0.1)

%% Clear out previous video frames
if strcmp(eType,'2 camera')
    delete([root filesep 'pike_tmp' filesep '*.TIF'])
    delete([root filesep 'marlin_tmp' filesep '*.TIF'])
else
    delete([root filesep 'video_tmp' filesep '*.TIF'])
end

%% Run the DAQ 

% Initiate sensor input channel(s)
if strcmp(eType,'camera + sensor')
    start(AI)
    trigger(AI)
end

% Step up the AO
putsample(AO,Vout)
disp('Camera triggered')

% Wait to stop running the input
if strcmp(eType,'camera + sensor')
    wait(AI,duration+1)
else
    pause(duration+1)
end

% Reset the AO 
putsample(AO,0)

% Grab data from memory
if strcmp(eType,'camera + sensor')
    data_in = getdata(AI);
else
    data_in = [];
end

% Delete DAQ objects for next run
if strcmp(eType,'camera + sensor')
    delete(AI)
end

delete(AO)

clear AI AO

%% Check video files

disp(' ')
disp('Giving Fire-i a moment to save the video . . .')

% Initialize 'tmp'
tmp = 0;

% Analyze filename format
if strcmp(eType,'2 camera')
    cIndx = strfind(capFileName,'$FI');
    
    % Check that camer index is there
    if isempty(cIndx) || cIndx==1
        error(['Capture file name needs to include $FI and it cannot be '...
               'the first item in the filename']);
    end 
    
    % Define filename prefixes
    prefix = [capFileName(1:cIndx-1)];
end

% Loop while video frames are writing to disk
while 1
    % Check only periodically
    pause(1)
    
    % List files in temporary directory
    if strcmp(eType,'2 camera')
        vFiles = dir([root filesep 'pike_tmp' filesep prefix '*.TIF']);
        vFiles1 = dir([root filesep 'marlin_tmp' filesep prefix '*.TIF']);
    else
        vFiles = dir([root filesep 'video_tmp' filesep '*.TIF']);
    end
    
    % Check progress
    if strcmp(eType,'2 camera') && (isempty(vFiles) && isempty(vFiles1))
        error(['Neither camera recorded -- You should have clicked on Start' ...
              ' in Fire-i before running this matlab program']);
        
    elseif strcmp(eType,'2 camera') && isempty(vFiles)
        error('Pike did not record')
        
    elseif strcmp(eType,'2 camera') && isempty(vFiles1)
        error('Marlin did not record')
        
    elseif isempty(vFiles) 
        error(['No video -- You should have clicked on Start' ...
              ' in Fire-i before running this matlab program']);
          
    elseif length(vFiles) >= numFrames
        break
        
    elseif length(vFiles) == tmp
        error(['You need to record at least ' num2str(numFrames+10) ...
           ' frames if your want to record at ' num2str(frame_rate) ...
           ' fps for ' num2str(duration) ' s.'])
    end
    
    % tmp is the number of files on the last loop
    tmp = length(vFiles);
    
    disp(['     . . . ' num2str(tmp) ' frames'])
end
disp(' Done !')

clear tmp

%% Prompt user to save data

answer = questdlg('Save recording?','','Yes','No','Yes');
if ~strcmp(answer,'Yes')
    return
end
  
% Generate directory name from individual, sequence numbers, & date
tmp1 = ['00' num2str(indivNum)];
tmp2 = ['00' num2str(seqNum)];
tmp3 = ['00' num2str(spirNum)];
dirName = ['R' tmp1(end-2:end) '-sp' tmp3(end-2:end) ...
           '-r' tmp2(end-2:end) '-' date];

% Try making the directory
[success,mess,messID] = mkdir([root filesep 'recordings' filesep dirName]);

% Browse to dir, if already there
if ~isempty(mess)
    but = questdlg('Overwrite existing data?',mess,'Yes','Create new',...
                   'Cancel','Yes');
               
    if strcmp(but,'Create new')
        cDir = uigetdir([root filesep 'recordings' filesep dirName],...
                            'Select where to save');
    elseif strcmp(but,'Yes')
        mess = [];
    else
        return
    end
end

% Otherwise make path for new dir
if isempty(mess)
    if strcmp(eType,'2 camera')
        cDir0 = [root filesep 'recordings' filesep dirName filesep 'pike'];
        cDir1 = [root filesep 'recordings' filesep dirName filesep 'marlin'];
        mkdir(cDir0)
        mkdir(cDir1)
    else
        cDir = [root filesep 'recordings' filesep dirName];
        mkdir(cDir)
    end
end

clear success mess messID tmp1 tmp2 tmp3

%% Move video files into recording directory

% Loop through each frame
for i = 1:numFrames
    if strcmp(eType,'2 camera')
        movefile([root filesep 'pike_tmp' filesep vFiles(i).name],cDir0);
        movefile([root filesep 'marlin_tmp' filesep vFiles1(i).name],cDir1);
    else
        movefile([root filesep 'video_tmp' filesep vFiles(i).name],cDir);
    end
end

% Delete remaining files
if strcmp(eType,'2 camera')
    delete([root filesep 'pike_tmp' filesep '*.TIF'])
    delete([root filesep 'marlin_tmp' filesep '*.TIF'])
else
    delete([root filesep 'video_tmp' filesep '*.TIF'])
end

%% Save recording data 

% Store in s structure
s.indiv_num   = indivNum;
s.spir_num    = spirNum;
s.seq_num     = seqNum;
s.exp_type    = eType;
s.ch_data     = data_in;
s.sample_rate = sample_rate_in;
s.frame_rate  = frame_rate;
s.pike_side   = camSide0;
s.marlin_side = camSide1;
s.file_names  = capFileName;

% Save s
if strcmp(eType,'2 camera')
    save([root filesep 'recordings' filesep dirName],'s')
else
    save([cDir filesep 'recording_data'],'s')
end

%% Prompt to run pix_acq

answer = questdlg('Run pixAcq on acquired video?','','Yes','No','Yes');

if ~strcmp(answer,'Yes')
    return
end

if strcmp(eType,'2 camera')
    pix_acq(vFiles(1).name,cDir0,frame_rate)
    pix_acq(vFiles(1).name,cDir1,frame_rate)
else
    pix_acq(vFiles(1).name,cDir,frame_rate)
end

clear answer

%% Prompt to run spir_analysis

answer = questdlg('Run spir_analysis on acquired video?','','Yes','No','Yes');

if ~strcmp(answer,'Yes')
    return
end

if strcmp(eType,'2 camera')
    spir_analysis([root filesep 'recordings' filesep dirName])
else
    spir_analysis([cDir filesep 'recording_data'])
end

clear answer