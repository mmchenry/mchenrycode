function compare_spir(pName)
% Compares two spiracle recordings 


%% Parameters

% Name of data files
fName = 'pixel_data.mat';

% High-pass cut-off frequency (Hz)
cut_high = .3;

% Tolerance for smoothing spline
tol = 1;

% Duration for box-car analysis (s)
anaDur = 10;


%% Get path of data file, load data

if nargin < 1
    pName = uigetdir(fName,'Select directory with pike and marlin');
end


%% Load pixel data

% From pike
load([pName filesep 'pike' filesep fName])

dP        = d.pixVal;
numFrames = d.numFrames;
frameRate = d.frameRate;

clear d

% From marlin
load([pName filesep 'marlin' filesep fName])

dM = d.pixVal;

% Check number of frames
if numFrames ~= d.numFrames
    error('Unequal number of frames between pike and marlin')
end

clear d

% Define time vector
t = [0:1/frameRate:(numFrames-1)/frameRate]';


%% Interactively select duration to analyze

a = dir([pName filesep 'interval_data.mat']);

if isempty(a)
    
    % High-pass filter pixel data
    dMf = butter_filt(dM,frameRate,cut_high,'high');
    dPf = butter_filt(dP,frameRate,cut_high,'high');
    
    % Normalize
    dMf = (dMf-mean(dMf))./std(dMf);
    dPf = (dPf-mean(dPf))./std(dPf);
    
    % Plot raw and normalized data --------------------
    figure;
    
    subplot(2,1,1)
    plot(t,dM,'r')
    hold on
    
    plot(t,dP,'b')
    legend('Marlin','Pike')
    
    xlabel('time (s)')
    ylabel('Raw mean pixel value')
    grid on
    
    subplot(2,1,2)
    plot(t,dMf,'r')
    hold on
    
    plot(t,dPf,'b')
    
    xlabel('time (s)')
    ylabel('Filtered/Normalized pixel value')
    grid on
    
    subplot(2,1,1)
    beep
    title('Select start point')
    
    % Prompt for interval question -----------------------
    ans = questdlg('Select interval to analyze?','',...
        'Yes','No, analyze all','Cancel','Yes');
    
    if strcmp(ans,'Cancel') || isempty(ans)
        return
        
    elseif strcmp(ans,'No, analyze all')
        tStart = 0;
        tEnd   = max(t);
        
    else
        
        % Select interval --------------------------------
        
        % Get point
        [x,y,b] = ginput(1);
        
        % Stop, if empty
        if isempty(x)
            return
        end
        
        tStart = x;
        
        % Plot selection
        hold on
        plot([tStart tStart],ylim,'k-')
        hold off
        
        clear x y b
        
        subplot(2,1,1)
        beep
        title('Select end point')
        
        % Get point
        [x,y,b] = ginput(1);
        
        % Stop, if empty
        if isempty(x)
            return
        end
        
        % Define ending time
        tEnd = x;
        
        % Plot selection
        hold on
        plot([tEnd tEnd],ylim,'k-')
        hold off
        
    end
    
    % Close window
    pause(.5)
    close
    
    clear x y b dMf dPf ans
    
    save([pName filesep 'interval_data'],'tStart','tEnd')
    
else
    disp(' ')
    disp('Loading previously-selected interval . . .')
    disp(' ')
    load([pName filesep 'interval_data'])
    
end

clear a
    

%% Trim, filter & normalize data

% Find indx for interval
idx = (t>=tStart) & (t<=tEnd);

% Trim data to interval
tf = t(idx);
tf = tf-tf(1);
dMf = dM(idx);
dPf = dP(idx);

% High-pass filter 
dMf = butter_filt(dMf,frameRate,cut_high,'high');
dPf = butter_filt(dPf,frameRate,cut_high,'high');

% Normalize
dMf = (dMf-mean(dMf))./std(dMf);
dPf = (dPf-mean(dPf))./std(dPf);

clear idx


%% Step through data for analysis

a = dir([pName filesep 'comparison_data.mat']);

if isempty(a)
    
    % Perform FFTs
    [f,r.pwrP] = fft_data(dPf,frameRate);
    [r.f,r.pwrM] = fft_data(dMf,frameRate);
    
    % Define index vector, initiate variables
    intrvl = floor(tf./anaDur)+1;
    i = 1;
    
    % Step through timeseries
    while 1
        
        % Break when completed all intervals
        if i > max(intrvl)
            break
        end
        
        % Pixel values for current interval
        cM = dMf(intrvl==i);
        cP = dPf(intrvl==i);
        cT = tf(intrvl==i);
        
        % Find delay (s)
        r.delay(i) = (length(cP) - find(xcorr(cM,cP)==max(xcorr(cM,cP))))/...
            frameRate;
        r.t_delay(i) = mean(cT);
        
        % Smoothing spline fit to data
        spM = spaps(cT,cM,tol);
        spP = spaps(cT,cP,tol);
        
        % Find zeros of splines (s)
        zM = fnzeros(spM)';
        zP = fnzeros(spP)';
        
        % Find mean period of openings (s)
        r.periodM(i) = 2*mean(diff(zM(:,1)));
        r.periodP(i) = 2*mean(diff(zP(:,1)));
        
        % Visualize correction from delay
        if 0
            subplot(3,1,1)
            plot(cT,cM,'r.',cT,cP,'b.');
            hold on
            fnplt(spM,'r')
            fnplt(spP,'b')
            grid on
            hold off
            
            subplot(3,1,2)
            plot(cT,cM,'r',cT-r.delay(i),cP,'b-');
            title('with delay correction:')
            grid on
            
            subplot(3,1,3)
            fnplt(fnder(spM,1),'r')
            hold on
            fnplt(fnder(spP,1),'b')
            hold off
            grid on
            
        end
        
        i = i+1;
        
        clear zM zP cM cP cTspM spP
    end
    
    % Save r
    save([pName filesep 'comparison_data'],'r')
    
    clear r i intrvl a
end


%% Display data

% Load 'r' vector
load([pName filesep 'comparison_data.mat'])

numPanels = 6;

figure
subplot(numPanels,1,1)
    plot(t-tStart,dM,'r',t-tStart,dP,'b');
    hold on
    
    yL = ylim;
    h = fill([tStart tEnd tEnd tStart]-tStart,...
             [yL(1) yL(1) yL(2) yL(2)],[.5 .5 .5]);
    set(h,'EdgeColor','none')
    set(h,'FaceAlpha',.3)
    clear h
    
    legend('Marlin','Pike')

    xlabel('time (s)')
    ylabel('Raw pix val')
    grid on

subplot(numPanels,1,2)
    plot(tf,dMf,'r')
    hold on

    plot(tf,dPf,'b')

    xlabel('time (s)')
    ylabel('Filtered pix val')
    grid on
    
subplot(numPanels,1,3)
    plot(r.t_delay,r.periodM,'ro-')
    hold on
    plot(r.t_delay,r.periodP,'bo-')
    hold off
    
    xlabel('time (s)')
    ylabel('Opening period (s)')
    grid on

subplot(numPanels,1,4)
    plot(r.t_delay,r.delay ,'ko-')
    
    xlabel('time (s)')
    ylabel('Pike delay behind marlin (s)')
    grid on
    
subplot(numPanels,1,5:6)
    plot(r.f,r.pwrM,'r',r.f,r.pwrP,'b')
    xlabel('Freq (Hz)')
    ylabel('Power')
    xlim([0 2])

    

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

function [f,P] = fft_data(y,sample_rate)
NFFT = 2^nextpow2(length(y)); % Next power of 2 from length of y
Y = fft(y,NFFT)/length(y);
f = sample_rate/2*linspace(0,1,NFFT/2+1);

P = 2*abs(Y(1:NFFT/2+1));

