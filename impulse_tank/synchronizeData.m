function [amp,daq,vid] = synchronizeData(mPath)
% Loads data from the daq and amplifier and synchronizes their timing such
% that t=0 occurs at the leading edge of the trigger signal

%%%%%%%%%%%%%%%%  SET PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%
% This constant used to convert the arbitrary amplifer position units into
% mm for continuity prediction
cal_const = 6.5274e-3; %mm /(amp position units)


%%%%%%%%%%%%%%%%%%%  SET PATHS  %%%%%%%%%%%%%%%%%%%%%%%%
motorFile = [mPath filesep 'motor.csv'];
daqFile   = [mPath filesep 'daqData.mat'];
videoFile = [mPath filesep 'videoData.mat'];


%% MOTOR KINEMATICS (from amplifier)

% The leading edge of the in9 signal is assumed to correspond to the
% global time = 0.  

%LOAD AND DEFINE =========================================
M       = importdata(motorFile,',');
in9     = M.data(:,4); %in V
t_amp   = M.data(:,1); %in s
pos_amp = M.data(:,2).*cal_const; %in mm

%CHECK ON COLUMN HEADINGS ================================
if ~strcmp(M.colheaders{1},'Time')
    error('First vector in the motor data should be: Time');
elseif ~strcmp(M.colheaders{2},'Actual motor position')
    error(['First vector in the motor data should be: ' ...
        'Actual motor position']);
elseif ~strcmp(M.colheaders{3},'Actual Motor Velocity')
    error(['First vector in the motor data should be: ' ...
        'Actual Motor Velocity']);
elseif ~strcmp(M.colheaders{4},'Input 9')
    error(['First vector in the motor data should be: ' ...
        'Input 9']);
end

%CALCULATE AMP_VEL =======================================
%filter data
cut_freq        = 100;
sample_rate     = 1./mean(diff(t_amp));
ff              = cut_freq./(sample_rate./2);
[B A]           = butter(2,ff,'low');
amp_pos_filt    = filtfilt(B,A,pos_amp); 
vel_amp         = diff(amp_pos_filt)./diff(t_amp);
t_vel_amp       = t_amp(2:end) - mean(diff(t_amp))/2;

clear cut_freq sample_rate B A amp_pos_filt


%CALCULATED VARIABLES ====================================
samplePeriod_amp    = mean(diff(t_amp));

%trigger onset time found as half sample period before rise above 0.05V:
trigStart_amp       = t_amp(find(in9>0.05,1)) - samplePeriod_amp/2;

%t_amp zeroed at time of intial rise at in9 on amp:
t_amp               = t_amp - trigStart_amp;
t_vel_amp           = t_vel_amp - trigStart_amp;

% Adjust position data to start at zero:
pos_amp = pos_amp - mean(pos_amp(1:50));

% Calculate start time of motor motion

threshAmp = 0.007; % Threshold position of motor when moving
threshDiff = 0.02; % Threshold change in position after which is ignored

%   1. estimate start of motion
iStart  = find(pos_amp>threshAmp,1)-1;
iEnd    = find(diff(pos_amp)>threshDiff,1);
% iEnd    = max(find((pos_amp(1:end-1)>threshAmp) & ...
%                 (diff(pos_amp)<threshDiff)));
    
%   2. find x-intercept for linear fit of position data 
c1 = polyfit(t_amp(iStart:iEnd),pos_amp(iStart:iEnd),1);
t_amp_mStart = -c1(2)/c1(1);


clear tHold iStart iEnd c1 threshAmp threshDiff

%% DAQ DATA  

% time = 0 for the daq is assumed to be the same instant as the global zero

%LOAD AND DEFINE ========================================
load(daqFile);

% Time vector from daq
t_daq       = daq.time; 

% Trigger signal from video
OUT1_daq    = daq.videoOUT1; 

% Trigger signal from arbitrary function generator
AFG_daq     = daq.AFG; 


%CALCULATED VARIABLES ====================================
sampleRate_daq  = 1./mean(diff(daq.time));

clear daq daqFile



%% VIDEO DATA  

% The time correction for video is calculated from the trailing edge of the 
% trigger signal (set to high value during fram exposures) and the total
% duration of the video

%LOAD AND DEFINE ========================================
load(videoFile);
rate_vid        = video.frameRate_fps;
dur_vid         = video.duration_frames;
postTrig_vid    = video.posTrig_frames;

%CALCULATED VARIABLES ====================================
% Use daq data and info about video trigger to calculate a
% time vector for the video that coincides with the daq data
exposurePeriod  = (1./rate_vid); 

if isempty(find(OUT1_daq<4.9)) || isempty(find(OUT1_daq>2))
    disp(' '); disp(['mPath = ' mPath]); disp(' ');
    warning(['No trigger signal in video -- setting zero time to half of ' ...
        'first frame after trigger']);
    
    t_frms  = (-(dur_vid-postTrig_vid):(postTrig_vid-1)) + 0.5;
    t_vid   = t_frms .* (1./rate_vid);
    clear frNum
else
    % video end time found from daq data (with high temporal resolution)
	% as a half sample prior to drop in OUT1 below 4.9V:
    endTime_vid = t_daq(min(find(OUT1_daq<4.9))) - 1./(2.*sampleRate_daq);

    % video start time calculated as (end time) - (video period):
    startTime_vid = endTime_vid - (dur_vid./rate_vid);

    % video time vector using frame rate and the start and end times:
    t_vid = startTime_vid:exposurePeriod:endTime_vid;

    % time vector corrected to instant in the middle of exposure
    t_vid = t_vid(1:end-1) + exposurePeriod/2; 

    % check that time vector length equal to number of frames in video
    if ~(length(t_vid)==dur_vid)
        error('Time vector not same length as video');
    end
end

clear endTime_vid startTime_vid video videoFile



%% STORE DATA  

daq.t      = t_daq;
daq.AFG    = AFG_daq;
daq.OUT1   = OUT1_daq;

amp.t       = t_amp;
amp.pos     = pos_amp;
amp.t_vel   = t_vel_amp;
amp.t_mStart= t_amp_mStart;
amp.vel     = vel_amp;
amp.in9     = in9;

vid.t       = t_vid;
vid.rate    = rate_vid;
vid.dur     = dur_vid;
vid.postTrg = postTrig_vid;



