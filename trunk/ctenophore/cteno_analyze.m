function cteno_analyze(dPath)


%% Parameters


%% Prompt for directory

% Prompt for first frame, if not given
if nargin < 1  
    dPath = uigetdir(pwd,'Choose directory with the data files');
    if isempty(dPath)
        return
    end
end

%% Load sequence data

% Look for existing sequence file
b = dir([dPath filesep 'seq_info.mat']);

% Check if present
if isempty(b)
    error('No sequence data file is present.');
end

% Load 'seq', structure of sequence data
load([dPath filesep 'seq_info']);

% Define parameters from seq
frame_rate = str2num(seq.frame_rate);
num_frames = seq.numFrames;

clear seq


%% Load plate data

% Look for existing plate_data file
b = dir([dPath filesep 'plate_data.mat']);

% Check if present
if isempty(b)
    error('No plate data exists for this sequence.');
end

% Load 'pl', structure of plate data
load([dPath filesep 'plate_data']);


%% Interactively select sequence for analysis

% Look for data file
b = dir([dPath filesep 'analysis_duration.mat']);
c = dir([dPath filesep 'beat_duration.mat']);

titles{1} = 'Select duration of analysis, then return';
titles{2} = 'Select duration of first full beat, then return';

if isempty(b) || isempty(c)
    
    for j = 1:2
        % Plot data
        if j == 1
            plot_combs(pl,num_frames,frame_rate)
            subplot(length(pl),1,1)
            
        end
        
        % Initialize variables
        xL = [nan nan];
        h  = [];
        
        % Add title
        beep
        title(titles{j})
        
        % Loop for collecting points
        while 1
            
            % Delete highlighted area
            if ~isempty(h)
                delete(h)
            end
            
            if j ==1
                
                % Plot area of x limits
                for i = 1:length(pl)
                    
                    subplot(length(pl),1,i)
                    
                    yL = ylim;
                    hold on
                    h(i) = fill([xL(1) xL(2) xL(2) xL(1)],...
                        [yL(1) yL(1) yL(2) yL(2)],'r');
                    alpha(h(i),0.3);
                    
                end
                
            else
%                 subplot(length(pl),1,1)
%                 
%                 for i = 1:length(pl)
%                     
%                     subplot(length(pl),1,i)
                    
                    yL = ylim;
                    xLold = frame_limits;
                    hold on
                    h(1) = fill([xLold(1) xLold(2) xLold(2) xLold(1)],...
                        [yL(1) yL(1) yL(2) yL(2)],'r');
                    alpha(h(1),0.3);
                    
                %end
                
                clear yL

                yL = ylim;
                %hold on
                h(2) = fill([xL(1) xL(2) xL(2) xL(1)],...
                    [yL(1) yL(1) yL(2) yL(2)],'b');
                alpha(h(2),0.3);
                %hold off
            end
            
            % Prompt for limit
            [x,y,b] = ginput(1);
            
            
            % If return
            if isempty(b)
                break
                
                % If left click
            elseif b == 1
                if isnan(xL(1))
                    xL(1) = x;
                else
                    xL(2) = x;
                end
                
                % If right click
            elseif b== 3
                if ~isnan(xL(2))
                    xL(2) = nan;
                else
                    xL = [nan nan];
                end
                
                %delete(h)
            end
        end
        
        if j == 1
            % Store results
            frame_limits = round([min(xL) max(xL)]);
            save([dPath filesep 'analysis_duration.mat'],'frame_limits');
        else
            
            if (min(xL) < frame_limits(1)) || (max(xL) > frame_limits(2))
                error('First beat must be within the duration selected')
            end
            
            beat_dur = round([min(xL) max(xL)]);
            save([dPath filesep 'beat_duration.mat'],'beat_dur');
        end
        
        close;
        
        if j == 1
            figure;
            
            % Define frame and time vectors
            frames = 1:num_frames;
            
            % Index for period to be analyzed
            idx = (frames>frame_limits(1)) & (frames<frame_limits(2));
            
            % Trim frames vector
            frames = frames(idx);
            
            angle = calc_angle(pl(1),idx);
            
            plot(frames,angle,'k-')
            h = [];
            
        end
        
    end
    
else
    % Load
    load([dPath filesep 'analysis_duration.mat'])
    load([dPath filesep 'beat_duration.mat'])
    
    disp(' ')
    disp('Loading from analysis_duration.mat & beat_duration.mat . . .')
    disp(' ')
end

clear b


%% 

figure

% Define frame and time vectors
frames = [1:num_frames]';

% Index for period to be analyzed
idx = (frames>frame_limits(1)) & (frames<frame_limits(2));

% Trim frames vector
frames = frames(idx);

% Loop through plates 
for i = 1 %:length(pl)
     
    % Calculate angle
    angle = calc_angle(pl(i),idx);
    
    % Low-pass filter
    angle_f = butter_filt(angle,frame_rate,15,'low'); 
    Dangle_f = diff(angle_f);
    
    % Extract first beat
    idx1    = (frames>beat_dur(1)) & (frames<beat_dur(2)); 
    frames1 = frames(idx1);
    angle1  = angle_f(idx1);
    tmp     = angle_f(find(idx1,1,'first')-1:find(idx1,1,'last'));
    Dangle1 = diff(tmp);
    tmp     = angle_f(find(idx1,1,'first')-1:find(idx1,1,'last'))+1;
    DDangle1= diff(diff(tmp));
    
    
    % First point 
    p1 = frames1(find(Dangle1==max(Dangle1),1,'first'));
    
    p2 = frames1(find((Dangle1<0) & (frames1>p1),1,'first'))-1;
    
    p3 = frames1(find((Dangle1==min(Dangle1)) & (frames1>p2),1,'first'))-1;
    
    % Step p1 back a half a power stroke cycle
    p1 = 2*p1-p2;
    
    % Step p3 back a half a recovery stroke cycle
    p3 = 2*p3-p2;
    
    subplot(3,1,1)
    plot(frames,angle,'k-',frames,angle_f,'r-')
    hold on
    plot([beat_dur(1) beat_dur(1)],ylim,'k--')
    plot([beat_dur(2) beat_dur(2)],ylim,'k--')
    plot([p1 p1],ylim,'b-')
    plot([p2 p2],ylim,'b--')
    plot([p3 p3],ylim,'g-')
    
    subplot(3,1,2)
    plot(frames(2:end)-0.5,Dangle_f,'r')
    hold on
    plot([p1 p1],ylim,'b-')
    plot([p2 p2],ylim,'b--')
    plot([p3 p3],ylim,'g-')
    grid on
    
    subplot(3,1,3)
    plot(frames(2:end-1),diff(Dangle_f),'r')
    hold on
    plot([p1 p1],ylim,'b-')
    plot([p2 p2],ylim,'b--')
    plot([p3 p3],ylim,'g-')
    grid on
end


% TODO: use these criteria to loop (1) down the combs plates and (2)
% through time on the first plate




function plot_combs(pl,num_frames,frame_rate)

% Set 'xlims' to unrealistic extremes
xlims = [10^5 0];

% Create figure window
figure

% Loop through plates
for i = 1:length(pl)
    
    % Define index where all coordinates exist
    idx = ~isnan(pl(i).tipY) & ~isnan(pl(i).baseY) & ~isnan(pl(i).angleY);
   
    % Calculate angle
    angle = calc_angle(pl(i),idx);
    
    % Define frame and time vectors
    frames = 1:num_frames;
    frames = frames(idx);
    %time   = (1/frame_rate).*frames(idx);
    
    % Plot
    subplot(length(pl),1,i)
    plot(frames, angle - min(angle),'k-')
    ylim([0 pi])
    
    % Update upper and lower limits to xaxis
    tmp = xlim;
    xlims = [min([xlims(1) tmp(1)]) max([xlims(2) tmp(2)])];
    
    clear tmp tip_vect base_vect angle frames
end

% Step though each plot to update x-axis
for i = 1:length(pl)
   subplot(length(pl),1,i)
   grid on
   xlim(xlims);
   ylabel(['angle (rad) of plate ' num2str(pl(i).plate_num)])
end

xlabel('frame number')


subplot(length(pl),1,1)

function angle = calc_angle(p,idx)
    % Angle of the tip wrt base in video frame coordinate system
    tip_vect = atan2((p.tipY(idx)-p.baseY(idx)),(p.tipX(idx)-p.baseX(idx)));
    
    % Angle of the body surface wrt base in video frame coordinate system
    base_vect = atan2((p.angleY(idx)-p.baseY(idx)),(p.angleX(idx)-p.baseX(idx)));
    
    % Difference btwn tip and base
    angle = tip_vect-base_vect;

    
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
