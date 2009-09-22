function d = headTracker
% Tracks the head of an adult fish that if swimming toward the left side of
% the video frame. vid_dir is the directory containing video tif files.
% use_centroid - (0 or 1) Use head centroid or posterior margin

% For kinematics acquisition figure:
vid_dir = '/Volumes/Docs/Projects/head_swimming/kinematic_data/Tales/ControlFish2_10';


use_centroid = 0;
action = 'acquire data'

%% General parameters

% Set path to m_file dir that contains this file
%m_path = '/Biology/Matlab/m_files';
m_path = '/Volumes/workgroup/m_file_library/kinematics/head_tracker';

% Specifies whether to display and frame analysis data during aquisition
visData   = 0;
visFrames = 1;


% Load default parameter values
load([m_path filesep 'default_params'])



%% Start case structure
switch action
    
    
%% Run start: initialize for a new video
case 'start'
    
% Prompt for video dir, save param
vid_dir = uigetdir(pwd,'Pick video directory');
if vid_dir==0
    return
else
    p.path = vid_dir;
end
save([m_path filesep 'default_params'],'p')

% Prompt whether to overwrite if parameters exist
if isfile('seq_params.mat',vid_dir)
      
    ButtonName = questdlg('Do you want to use previous sequence parameters?', ...
        'Data file found', ...
        'Yes', 'No, reset parameters', 'Yes');
    close
    
    switch ButtonName,
        case 'No, reset parameters',
            headTracker('set parameters',vid_dir)
            headTracker('acquire data',vid_dir)
        case 'Yes',
            headTracker('acquire data',vid_dir)
    end % switch
    
else %if no previous data
    headTracker('set parameters',vid_dir)
    headTracker('acquire data',vid_dir)
end %isfile

   
%% Define sequence parameter values
case 'set parameters'
    

    
    p.path = vid_dir;
    clear vid_dir
    
    % Prompt for im, if file absent
    if ~isfile(p.fname,p.path)
        [p.fname, pathName, filterIndex] = uigetfile([p.path filesep '*.tif'], ...
            'Select video frame in this directory');
        if filterIndex==0
            return
        end

        pathName = pathName(1:end-1);
        
        if ~strcmp(pathName,p.path)
            error('You cannot change directories')
        end
    end
    
    % Load im
    im = imread([p.path filesep p.fname]);
    
    % Find the threshold value that works for a new movie
    waitfor(threshFinder(im,p))
    load([p.path filesep 'seq_params.mat'])

    
    % Find a calibration constant
    [fName_cal, pathName_cal, filterIndex] = uigetfile('*.tif', ...
        'Select calibration image');
    if filterIndex==0
        return
    end
    im_cal = imread([pathName_cal filesep fName_cal]);
    p.calconst = calibrate(im_cal);
    save([p.path filesep 'seq_params.mat'],'p')
    

    % Find the length of the fish's head (pix)
    [x,y] = choosePoints(im,2,1,'Choose head length');
    p.headlength = ((x(2)-x(1))^2 + (y(2)-y(1))^2)^0.5;
    save([p.path filesep 'seq_params.mat'],'p')
    
    
    % Select region of interest (roi)
    figure;
    imshow(im)
    
    ButtonName = questdlg('Do you want to select a region of interest?', ...
        'ROI?', ...
        'Yes', 'No, use whole frame', 'Yes');
    close
    
    switch ButtonName,
        case 'No, use whole frame',           
            p.col_min   = 1;
            p.col_max   = size(im,2);
            p.row_min   = 1;
            p.row_max   = size(im,1);

        case 'Yes',
            [xR,yR] = choosePoints(im,2,2,'Choose ROI');            
            p.col_min   = round(min(xR));
            p.col_max   = round(max(xR));
            p.row_min   = round(min(yR));
            p.row_max   = round(max(yR));         
    end % switch
    
    
    % Get other parameters
    files = giveTiffStack(p.path,p.fname);
    prompt={'Frame rate (per sec)','Linear units','Start frame','End frame',...
        'Frame increment','Flow tank motor setting'};
    name='Other parameters';
    numlines=1;
    defaultanswer={num2str(p.framerate),'cm','1',num2str(length(files)),'1',...
                    '20'};
 
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    
    p.framerate  = str2num(answer{1});
    p.units      = answer{2};
    p.startframe = str2num(answer{3});
    p.endframe   = str2num(answer{4});
    p.framestep  = str2num(answer{5});
    p.motorset   = str2num(answer{6});
    
    % Store p
    save([p.path filesep 'seq_params.mat'],'p')
    save([m_path filesep 'default_params.mat'],'p')
    
    
%% Acquire coordinate data from each frame
case 'acquire data'
    
% Prompt for dir, if not given
    
    load([vid_dir filesep 'seq_params.mat'])
    
%         startFrame = p.startframe;
%         endFrame   = p.endframe;
%         frameStep  = p.framestep;
    %files = giveTiffStack(vid_dir,p.fname);
    files = dir([vid_dir filesep '*.tif']);
    j     = 1;
    
%     % Determine the index iNum of the ending of the file name that
%     % includes the frame number
%     name = files(1).name;
%     iZeros = find(name=='0');
%     if max(diff(find(name=='0'))) > 1
%         firstZero = (iZeros(max(find(diff(iZeros)>1))+1));
%     else
%         firstZero = min(find(name=='0'));
%     end
%     iNum        = firstZero:length(name);
    
    % define start of file name
    nameHead    = files(1).name(1:20);
    
       
    % Define roi
    roi.c = p.col_min:p.col_max;
    roi.r = p.row_min:p.row_max;      
        
    if visFrames
        figure
        set(gcf,'DoubleBuffer','on')
    end
    
    % Loop through each image in vid_dir
    for i = 397 %1:length(files)
        
        frNum = str2num(files(i).name(length(nameHead)+1:end-4));
        
        % Read and crop image
        im = imread([vid_dir filesep files(i).name]);
        im = im(roi.r,roi.c);
        
        % Run 'get_head'
        [xTip,yTip,xCtr,yCtr]  = get_head(im,p.headlength,p.tVal,visData,...
                                          use_centroid);
        % Visualize                              
        if visData
            subplot(2,2,1);
            title(['Frame ' num2str(i) ' out of ' num2str(length(files))]);
            subplot(2,2,4)
        end
        
        %Visualize tracking
        if visFrames
            warning off
            imshow(im)
            hold on
            plot([xTip xCtr],[yTip yCtr],'r-',...
                xTip,yTip,'or')
            title(['Frame number ' num2str(frNum)]);
            hold off
            pause(.01);
            warning on
        end
        
        %Store data
        
        d.framenum(j,1) = frNum;
        d.filename{j} = files(i).name;
        d.xTip(j,1) = xTip;
        d.yTip(j,1) = yTip;
        d.xCtr(j,1) = xCtr;
        d.yCtr(j,1) = yCtr;
        d.use_centroid(j,1) = use_centroid;
        
        %disp(['Frame ' num2str(i) '/' num2str(length(files))]);
        j = j+1;
        
    end
    
    if use_centroid
        save([vid_dir filesep 'coord_data'],'d');   
    else
        save([vid_dir filesep 'coord_data_posterior'],'d');   
    end
    
end %case



%% FUNCTIONS

function yes = isfile(fname,dirname)
% Checks if file is present in given directory
a = dir(dirname);
yes = 0;
for i = 3:length(a)
    if strcmp(a(i).name,fname)
        yes = 1;
        return
    end
end


function data_filtered = butterworth(data,sample_rate,cutfreq,type)
%Returns data filtered by a butterworth filter
%  data - vector of data for filtering
%  sample_rate - rate of sampling for data (must have equivalent intervals)
%  cutfreq - cut off frequency (must be 2 values for bandstop)
%  type - 'low' for lowpass (default), 'high' for highpass, 'stop' for bandstop

if nargin < 4
    type = 'low';
end

if strcmp(type,'stop') && ~(length(cutfreq)==2)
    error('Stop pass must have a two-element cutfreq')
end 

ff              = cutfreq./(sample_rate./2);
[B A]           = butter(2,ff,type);
data_filtered   = filtfilt(B,A,data);     % Filtered data


function files = giveTiffStack(mpath,fname)
% Returns a structure with info on the tiff stack

[pathstr,name,ext,versn]    = fileparts(fname);

clear pathstr ext versn

% % Determine the index iNum of the ending of the file name that 
% % includes the frame number
% iZeros = find(name=='0');
% if max(diff(find(name=='0'))) > 1
%     firstZero = (iZeros(max(find(diff(iZeros)>1))+1));
% else
%     firstZero = min(find(name=='0'));
% end
% iNum        = firstZero:length(name);

% define start of file name
nameHead    = name(1:20);


% set up for loop
a           = dir([mpath filesep '*.tif']);
startNum    = str2num(a(1).name(length(nameHead)+1:end-4));
tNum        = startNum;
j           = 1;
while 1==1
    nameEnd     = [num2str(zeros(1,...
                           length(name)-length(nameHead)-...
                           length(num2str(tNum)))) num2str(tNum)];
    tempName    = [nameHead nameEnd(find(~(nameEnd==' ')))];
    isFile      = 0;
    %Step through file list to see if file exists:
    for i = (tNum-startNum)+1:length(a)
        [pathstr,name,ext,versn]    = fileparts(a(i).name);
        if ~(min(ext=='.')) & (length(name)>=length(tempName))
            if strcmp(name(1:length(tempName)),tempName)
                files(j).name   = tempName;
                j               = j + 1;
                tNum            = tNum + 1;
                isFile          = 1;
                break
            end
        end
    end
    if ~isFile, break; end
end
if j==1, files = []; end