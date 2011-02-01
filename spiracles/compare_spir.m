function compare_spir(pName)
% Compares two spiracle recordings 


%% Default parameters

% Name of data files
fName = 'pixel_data.mat';

% High-pass cut-off frequency (Hz)
p.cut_high = .3;

% Tolerance for smoothing spline
p.tolM = 1;
p.tolP = 1;

% Duration for box-car analysis (s)
p.anaDur = 15;

% Step through analysis
stepThrough = 1;

% Whether to invert pixel intensity values
p.invert = 'no';


%% Get path of data file, load data

if nargin < 1
    pName = uigetdir(fName,'Select directory with pike and marlin');
    if pName==0
        return
    end
end

%% Prompt about running

a  = dir([pName filesep 'interval_data.mat']);
a2 = dir([pName filesep 'include_in_pool.mat']);
a3 = dir([pName filesep 'exclude_from_pool.mat']);

if ~isempty(a)
    
    if ~isempty(a2)
        disp(' ')
        disp(['Sequence approved (file include_in_pool.mat exists), '...
            'visualizing results . . .'])
        
    elseif ~isempty(a2)
        disp(['Sequence rejected (file exclude_from_pool.mat exists), '...
            'visualizing results . . .'])
        
    else
        
        but = questdlg('Data exist','Warning!',...
            'Overwrite all','Keep interval',...
            'Just visualize','Overwrite all');
        
        if isempty(but)
            return
            
        elseif strcmp(but,'Overwrite all')
            delete([pName filesep 'comparison_data.mat'])
            delete([pName filesep 'interval_data.mat'])
            
            a = dir([pName filesep 'ana_params.mat']);
            if ~isempty(a)
                delete([pName filesep 'ana_params.mat'])
            end
            
            a = dir([pName filesep 'comparison_data.mat']);
            if ~isempty(a)
                delete([pName filesep 'comparison_data.mat'])
            end
            
        elseif strcmp(but,'Keep interval')
            delete([pName filesep 'comparison_data'])
            delete([pName filesep 'ana_params.mat'])
            
        elseif strcmp(a,'Just visualize')
            % Do nothing
            
        end
        
    end
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
    dMf = butter_filt(dM,frameRate,p.cut_high,'high');
    dPf = butter_filt(dP,frameRate,p.cut_high,'high');
    
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

% Load parameter file, if one exists
a2 = dir([pName filesep 'ana_params.mat']);

if ~isempty(a2)
    % Load 'p' structure
    load([pName filesep 'ana_params.mat'])
    
end

% Find indx for interval
idx = (t>=tStart) & (t<=tEnd);

% Trim data to interval
tf = t(idx);
tf = tf-tf(1);
dMf = dM(idx);
dPf = dP(idx);

p.cut_high=.1;

% High-pass filter 
dMf = butter_filt(dMf,frameRate,p.cut_high,'high');
dPf = butter_filt(dPf,frameRate,p.cut_high,'high');

if 0
    figure;
    subplot(2,1,1)
        h = plot(tf,dM(idx)-mean(dM(idx)),'r',tf,dMf,'r');
        set(h(2),'LineWidth',1.5)
    subplot(2,1,2)
        plot(tf,dP(idx)-mean(dP(idx)),'b',tf,dPf,'b');
        set(h(2),'LineWidth',1.5)
end

% Normalize filtered
dMf = (dMf-mean(dMf))./std(dMf);
dPf = (dPf-mean(dPf))./std(dPf);


%% Step through data for analysis

a = dir([pName filesep 'comparison_data.mat']);

if isempty(a)
    
    % Perform FFTs
    [f,r.pwrP] = fft_data(dPf,frameRate);
    [r.f,r.pwrM] = fft_data(dMf,frameRate);
    
    % Load parameter file, if one exists
    a2 = dir([pName filesep 'ana_params.mat']);
    if isempty(a2)
        stepThrough = 1;
        
    else
        % Load 'p' structure
        load([pName filesep 'ana_params.mat'])
        
    end
    
    if stepThrough
        figure
    end
    
    % Define index vector, initiate variables
    i = 1;
    dspM = [];
    dspP = [];
    
    % Step through timeseries
    while 1
        
        % Loop that checks interval duration
        while 1
            % Define interval
            intrvl = floor(tf./p.anaDur)+1;
            
            % Break when completed all intervals
            if i > max(intrvl)-1
                break
            end
            
            % Pixel values for current interval
            cM = dMf(intrvl==i);
            cP = dPf(intrvl==i);
            cT = tf(intrvl==i);
            
            % Invert, if requested
            if strcmp(p.invert,'pike')
                cP = -cP;
            elseif strcmp(p.invert,'marlin')
                cM = -cM;
            elseif strcmp(p.invert,'both')
                cM = -cM;
                cP = -cP;
            end
                
            
            % Smoothing spline fit to data
            spM = spaps(cT,cM,p.tolM);
            spP = spaps(cT,cP,p.tolP);
            
            % Find zeros of splines (s)
            zM = fnzeros(spM)';
            zP = fnzeros(spP)';
            
            % Choose only zeros on leading edge of wave
            zMn = [];
            zPn = [];
            
            for j = 1:length(zM)
                if fnval(fnder(spM),zM(j))>0
                    zMn = [zMn; zM(j)];
                end
            end
            
            for j = 1:length(zP)
                if fnval(fnder(spP),zP(j))>0
                    zPn = [zPn; zP(j)];
                end
            end
            
            zM = zMn;
            zP = zPn;
            
            clear zMn zPn
            
            % Check result
            if  length(zM)>1 && length(zP)>1 
                break
                
            else
                disp(' ')
                beep
                disp(['Fewer than 2 peaks in interval! Increasing '...
                      'duration from ' num2str(p.anaDur) 's to ' ...
                      num2str(1.5*p.anaDur) 's'])
                disp(' ')
                
                p.anaDur = 1.5*p.anaDur;
                
                if p.anaDur > tf
                    error('Interval exceeds total duration')
                end
            end

        end
        
        % Break when completed all intervals
        if i > max(intrvl)-1
            break
        end
        
        % Find delay (s)
        r.delay(i) = (length(cP) - find(xcorr(cM,cP)==max(xcorr(cM,cP))))/...
            frameRate;
        r.t_delay(i) = mean(cT);
            
        % Find mean period of openings (s)
        r.periodM(i) = mean(diff(zM(:,1)));
        r.periodP(i) = mean(diff(zP(:,1)));
        
        % Evaluate splines for comparison with data
        dspM = [dspM; fnval(spM,cT)];
        dspP = [dspP; fnval(spP,cT)];
        
        % Set constants for inversion
        if strcmp(p.invert,'marlin')
            M_in = -1;
            P_in = 1;
        elseif strcmp(p.invert,'pike')
            P_in = -1;
            M_in = 1;
        elseif strcmp(p.invert,'both')
            P_in = -1;
            M_in = -1;
        else
            M_in = 1;
            P_in = 1;
        end
        
        
        % Visualize
        if stepThrough
            
            subplot(7,1,1)
            plot(tf,M_in.*dMf,'r',tf,P_in.*dPf,'b');
            hold on
            yL = ylim;
            h = fill([min(cT) max(cT) max(cT) min(cT)],...
                [yL(1) yL(1) yL(2) yL(2)],[.5 .5 .5]);
            set(h,'EdgeColor','none')
            set(h,'FaceAlpha',.3)
            clear h
            hold off
            
            subplot(7,1,2:3)
            h = plot(cT,cM,'r-');
            hold on
            h = plot(cT,fnval(spM,cT),'r-');
            set(h,'LineWidth',1.5)
            for j = 1:length(zM)
                plot([zM(j) zM(j)],ylim,'k-')
            end
            plot(xlim,[0 0],'k--')
            hold off
            title('Interval')
            xlabel('time (s)')
            ylabel('Normalized pixel intensity')
            legend('Marlin data','spline fit','zero cross')
            
            subplot(7,1,4:5)
            h = plot(cT,cP,'b-');
            hold on
            h = plot(cT,fnval(spP,cT),'b-');
            set(h,'LineWidth',1.5)
            for j = 1:length(zP)
                plot([zP(j) zP(j)],ylim,'k-')
            end
            plot(xlim,[0 0],'k--')
            hold off
            title('Interval')
            xlabel('time (s)')
            ylabel('Normalized pixel intensity')
            legend('Pike Data','spline fit','zero cross')
            
            subplot(7,1,6:7)
            plot(cT,cM,'r',cT-r.delay(i),cP,'b-');
            xlabel('time (s)')
            ylabel('Normalized pixel intensity')
            title('Measurements with delay correction:')
            grid on
            
            but = questdlg('How do the data look ?','',...
                     'Advance interval','Change parameter values',...
                     'Great, skip to end','Advance interval');
                 
            if isempty(but)
                return
                
            elseif strcmp(but,'Change parameter values')
                
                % Prompt for parameter values
                prompt = {'Interrogation interval (s)',...
                    'Marlin smoothing tolerance (larger = smoother)',...
                    'Pike smoothing tolerance (larger = smoother)',...
                    'Invert? ("no","pike","marlin","both")'};
                
                answer = inputdlg(prompt,'Parameters',1,...
                         {num2str(p.anaDur),...
                          num2str(p.tolM),num2str(p.tolP),p.invert});
                
                p.anaDur = str2num(answer{1});
                p.tolM = str2num(answer{2});
                p.tolP = str2num(answer{3});
                
                % Check/set invert input
                if ~strcmp(answer{4},'no') && ~strcmp(answer{4},'pike')...
                        && ~strcmp(answer{4},'marlin')...
                        && ~strcmp(answer{4},'both')
                    warning('invalid entry for "invert", setting to "no"')
                    p.invert = 'no';
                else
                    p.invert = answer{4};
                end
                
                % Rerun current frame
                i = i-1;
                
                clear prompt answer
            
            elseif strcmp(but,'Great, skip to end') 
                stepThrough = 0;
                
            end
            
            clear but 
        end
        
        i = i+1;
        
        clear zM zP cM cP cTspM spP
    end 
    
    % Save r & p
    save([pName filesep 'comparison_data'],'r')
    save([pName filesep 'ana_params.mat'],'p')
    
    clear r i intrvl a p
end


%% Display data 

% Load 'r' vector
load([pName filesep 'comparison_data.mat'])

% Load 'p' vector
load([pName filesep 'ana_params.mat'])

numPanels = 6;

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

figure
subplot(numPanels,1,1)
    plot(t-tStart,M_in.*dM,'r',t-tStart,P_in.*dP,'b');
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
    plot(tf,M_in.*dMf,'r')
    hold on

    plot(tf,P_in.*dPf,'b')

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

    
%% Prompt for revision

a2 = dir([pName filesep 'include_in_pool.mat']);
a3 = dir([pName filesep 'exclude_from_pool.mat']);

if isempty(a2) && isempty(a3)
    
    but = questdlg('How do the data look ?','',...
        'Great!','Re-run analysis',...
        'Poor -- do not include','Great!');
    
    if isempty(but)
        return
        
    elseif strcmp(but,'Great!')
        include_in_pool = 1;
        save([pName filesep 'include_in_pool.mat'],'include_in_pool');
        
        but2 = questdlg('Are the periods the same for the two spiracles',...
                        '','Yes','No','Yes');
        if isempty(but2)
            return
            
        elseif strcmp(but2,'No')
            save([pName filesep 'unequal_periods.mat'],'include_in_pool');
            
        else 
            save([pName filesep 'equal_periods.mat'],'include_in_pool');
            
        end
        
        clear but2
            
        
    elseif strcmp(but,'Re-run analysis')
        close
        %     delete([pName filesep 'comparison_data'])
        %     delete([pName filesep 'ana_params.mat'])
        %     delete([pName filesep 'interval_data.mat'])
        compare_spir(pName);
        
    elseif strcmp(but,'Poor -- do not include')
        include_in_pool = 0;
        save([pName filesep 'exclude_from_pool.mat'],'include_in_pool');
        
    end
    
    clear but
    
end





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

