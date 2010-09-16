function spir_animate(pName)
% Animates data collected with pixAcq
%
% pName - path to images and data file
%
% Note: requires that the data file created by pixAcq be in the same
% directory as the image files

%% Parameters

% Name of data file
fName = 'pixel_data.mat';

% High-pass cut-off frequency (Hz)
cut_high = 0.3;

% Rate to play video
%play_rate = 20;

frame_skip = 1;

winSize = 750;

% Proportion of sequence to create movie:
duration = 0.25;

% Period over which to visualize the signal (s)
prd = 30;


%% Get path of data file, load data

if nargin < 1
    pName = uigetdir(fName,'Select directory with pike and marlin video');
    if pName==0
        return
    end
end

% Load 'p' vector
load([pName filesep 'ana_params.mat'])

% Set constants for inversion
if strcmp(p.invert,'marlin')
    M_in = -1;
    P_in = 1;   
elseif strcmp(p.invert,'pike')
    P_in = -1;
    M_in = 1;    
else
    M_in = 1;
    P_in = 1;
end

% Load pike data
load([pName filesep 'pike' filesep fName])
aP = dir([pName filesep 'pike' filesep d.fileName{1}(1:2) '*']);
dP = d;
dP.pixVal = P_in.*dP.pixVal;
clear d

% Load marlin data
load([pName filesep 'marlin' filesep fName])
aM = dir([pName filesep 'marlin' filesep d.fileName{1}(1:2) '*']);
dM = d;
dM.pixVal = M_in.*dM.pixVal;

frameRate = d.frameRate;
t = d.frameNum./d.frameRate;

clear d p M_in P_in


%% Make figure

f = figure;
set(f,'DoubleBuffer','on')
%set(f,'WindowStyle','docked')
set(f,'WindowStyle','normal')
%set(f,'WindowStyle','modal')

%% Process pixel data

% For pike
pixP = butter_filt(dP.pixVal,dP.frameRate,cut_high,'high');
pixP = (pixP-mean(pixP))./std(pixP);

pixM = butter_filt(dM.pixVal,dM.frameRate,cut_high,'high');
pixM = (pixM-mean(pixM))./std(pixM);


%% Plot data

subplot(3,2,5:6)
plot(t,pixM,'r',t,pixP,'b')
xlabel('Time (s)')
ylabel('Mean pixel intensity')
hold on
yVals = ylim;


%% Loop through frames

% Initialize movie file
if ispc
    aviobj = avifile([pName filesep 'data_movie'],'compression','Cinepak');
else
    aviobj = avifile([pName filesep 'data_movie'],'compression','None');
end

for i = 1:frame_skip:duration*length(dM.frameNum)
    
    tic
    
    % Define current time
    tc = dM.frameNum(i)./dM.frameRate;
    
    % Marlin: Read current frame, roi
    IM = imread([pName filesep 'marlin' filesep aM(i).name]);
    
    roiIM = roipoly(IM,dM.xROI{i},dM.yROI{i});
    bIM = uint8(255.*ones(size(IM,1),size(IM,2)));
    bIM(roiIM(:))=IM(roiIM(:));
    
    % Pike: Read current frame, roi
    IP = imread([pName filesep 'pike' filesep aP(i).name]);
    roiIP = roipoly(IP,dP.xROI{i},dP.yROI{i});
    bIP = uint8(255.*ones(size(IP,1),size(IP,2)));
    bIP(roiIP(:))=IP(roiIP(:));
      
    % Mark current time on graph
    figure(f)
    set(f,'Position',[1 winSize winSize winSize])
    subplot(3,2,5:6)
    h1 = plot([tc tc],yVals,'k-');
    xlim([tc-prd/4 tc+3*prd/4])
    
    % Display current marlin video frame
    subplot(3,2,[1,3])
    warning off
    %imshow(bIM);
    imshow(bIM)
    axis([min(dM.xROI{i}) max(dM.xROI{i}) ...
          min(dM.yROI{i}) max(dM.yROI{i})])
        
    %hold on
    warning on
    %plot(dM.roiX,dM.roiY,'r-')
    %plot(dM.xROI{i},dM.yROI{i},'r-')
    %hold off
    h = title('Marlin');
    set(h,'Color','r')
    
    % Display current pike video frame
    subplot(3,2,[2,4])
    warning off
    %imshow(bIP);
    imshow(bIP)
    %hold on
    warning on
    %plot(dP.roiX,dP.roiY,'b-')
    %plot(dP.xROI{i},dP.yROI{i},'b-')
    %hold off
    h = title('Pike');
    set(h,'Color','b')
    axis([min(dP.xROI{i}) max(dP.xROI{i}) ...
          min(dP.yROI{i}) max(dP.yROI{i})])
    
    %         % Pause, then proceed
    %         elTime = toc;
    %         delay = max([0 (1/play_rate)-elTime]);
    %         pause(delay)
    
    Fr = getframe(f);
    aviobj = addframe(aviobj,Fr);
    
    % Delete line
    delete(h1)
    
    % Clear variables
    clear t eTime delay Fr
    
end

close(f);
aviobj = close(aviobj);


function data_filtered = butter_filt(data,sample_rate,cut_freq,type) 
% High-pass or low-pass butterworth filter

% All frequency values are in Hz.

% Nyquist freq.
Nqst = sample_rate/2;   

% Calculate stopband frequency
if strcmp(type,'high')
    stop_freq = max([(cut_freq - Nqst/10) .01]);  

elseif strcmp(type,'low')
    stop_freq = min([(cut_freq + Nqst/10) (Nqst-.01)]); 
 
end

% Stopband Attenuation (dB)
Astop = 10;

% Passband Ripple (dB)
Apass = 1;   

% Normalise the cutoff freq. to the Nyquist freq
Wp = cut_freq/Nqst;

% Normalise the stoppass freq. to the Nyquist freq
Ws = stop_freq/Nqst;

% Check cutoff
if (Wp > 1) || (Ws > 1)
    error('Cutoff or bandpass frequencies cannot exceed the Nyquist freq')
end

% Calculate the order from the parameters using BUTTORD.
[N,Fc] = buttord(Wp, Ws, Apass, Astop);    
    
% Calculate the zpk values using the BUTTER function.
[B A] = butter(N, Fc, type);

% Plot frequency reponse
%freqz(B,A,512,sample_rate); 

% Filter the data
data_filtered   = filtfilt(B,A,data); 

