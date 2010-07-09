function angle_acquire(imPath,f_name)


%% Parameters

num_digits = 6;     % Digits at end of image filenames
name_prefix = 'cteno';


%% Define directories

% Prompt for first frame, if not given
if nargin < 2  
    [f_name,imPath,fIndex] = uigetfile({'*.tif','*.jpg'},...
        'Choose first image in sequence');
    if ~fIndex
        return
    end
end

%% Acquire sequence information

% Look for existing seq_info file and load file
a = dir([imPath filesep 'seq_info.mat']);

if isempty(a)
    %output error message
    error('No sequence data exists for this sequence.');
else
load([imPath filesep 'seq_info']);
end

% Look for existing plate_data file and load file
b = dir([imPath filesep 'plate_data.mat']);

if isempty(b)
    %output error message
    error('No plate data exists for this sequence.');
else
    load([imPath filesep 'plate_data']);
end

xlims = [10^5 0];
for i=1:length(pl)
    
    idx = ~isnan(pl(i).tipY) & ~isnan(pl(i).baseY) & ~isnan(pl(i).angleY);
   
    tip_vect = atan2((pl(i).tipY(idx)-pl(i).baseY(idx)),(pl(i).tipX(idx)-pl(i).baseX(idx)));
    base_vect = atan2((pl(i).angleY(idx)-pl(i).baseY(idx)),(pl(i).angleX(idx)-pl(i).baseX(idx)));
    
    angle = tip_vect-base_vect;
    
    numFrames = 1:seq.numFrames;
    time = (1/125).*numFrames(idx);
    
    subplot(length(pl),1,i)
    plot(time, angle - min(angle))
    ylim([0 3])
    grid on
    tmp = xlim;
    xlims = [min([xlims(1) tmp(1)]) max([xlims(2) tmp(2)])];
    
    clear time 
    %pause
end

for i = 1:length(pl)
   subplot(length(pl),1,i)
   %xlim(xlims);
   xlim([0 4.5])
end

%figures
% figure
% subplot(5,1,1)
% plot(time(1),angle(1))
% ylabel('Plate 1')
% subplot(5,1,2)
% plot(time(2),angle(2))
% ylabel('Plate 2')
% subplot(5,1,3)
% plot(time(3),angle(3))
% ylabel('Plate 3')
% subplot(5,1,4)
% plot(time(4),angle(4))
% ylabel('Plate 4')
% subplot(5,1,5)
% plot(time(4),angle(5))
% ylabel('Plate 5')
% xlabel('time (s)')
