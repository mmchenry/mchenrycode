function acquire_strain(vPath)
% Analyzes the video from a cupula material testing experiment

%% Parameters

% Whether to enhance contrast
enhance_contrast = 1;

% Number of digits for frame numbers
num_dig = 6;

% num_pix is the number of pixels averaged over to find an edge
num_pix = 10;

% window size to assess position of pipette tip
pix_window = 20;

% Default values for short and long sequences
fps_short = 1000;
fps_long  = 125;

% Interval for skipping thru frames in data preview
fSkip = 10;

% Default smoothing tolerance
smooth_def = 1;


%% Find images

% Prompt for directory
if nargin < 1
   [vPath] = uigetdir(pwd,'Choose video directory');
end

% Check for data file
if ~isempty(dir([vPath filesep 'edge_data.mat'])) || ...
   ~isempty(dir([vPath filesep 'seq_data.mat']))

    but = questdlg('Sequence already started?', ...
                         'Warning!', 'Start over','Use seq_data', 'Stop', 'Stop');
    if strcmp(but,'Stop')
        return
        
    elseif strcmp(but,'Start over')
        delete([vPath filesep 'seq_data.mat'])
        
        
    end
end

% List tif files
aTMP = dir([vPath filesep '*.tif']);

% Check input
if isempty(aTMP)
    error('No tif files in the selected directory')
end

% Define file prefix
fname_pre = aTMP(1).name(1:end-4-num_dig);

% List of all tiff files with the prefix
p.files = dir([vPath filesep fname_pre '*.tif']);

clear aTMP


%% Prompt for sequence info

if isempty(dir([vPath filesep 'seq_data.mat']))
    % Get info on experiment
    but = questdlg('What type of experiment?','','Short','Long','Short');
    
    
    % Prompt for seq parameters
    if isempty(but)
        
        return
        
    elseif strcmp(but,'Short')
        
        p.exp_type = 'short';
        answer = inputdlg({...
            'Frame number to start: ',...
            'Frame number to end: ',...
            'Frame rate (fps):',...
            'Mag cube:',...
            'Temperature (deg C):'...
            'Inner diameter (microns)',...
            },'',1,...
            {'1',num2str(length(p.files)),num2str(fps_short),'1.6',...
             '20','20'});
        
    elseif strcmp(but,'Long')
        
        p.exp_type = 'long';
        answer = inputdlg({...
            'Frame number to start:',...
            'Frame number to end:',...
            'Frame rate (fps):',...
            'Mag cube:',...
            'Temperature (deg C):'...
            'Inner diameter (microns)',...
            },'',1,...
            {'1',num2str(length(p.files)),num2str(fps_long),'1.6',...
             '20','20'});
        
    end
    
    % Check if cancelled
    if isempty(answer)
        return
    end
    
    p.start_frame  = str2num(answer{1}); 
    p.end_frame    = str2num(answer{2});
    p.frame_rate   = str2num(answer{3});
    p.temp         = str2num(answer{5});
    p.inner_diam   = str2num(answer{6}) * 10^-6;
    
    % Check mag cube input
    if strcmp(answer{4},'1.6')
        p.magcube = 1.6;
        
    elseif strcmp(answer{4},'2.5')
        p.magcube = 2.5;
        
    elseif strcmp(answer{4},'1.0') || strcmp(answer{4},'1') 
         p.magcube = 1.0;
         
    else
        error('magcube entry must be "1.6", "2.5", or "1')
        
    end
        
    % p.end_frame = 400;
    % p.frame_rate = 1000;
    
    clear fps_short fps_long fend_short fend_long
    
    
%% Select tracking line
    
    if strcmp(p.exp_type,'short')
        cFrame = p.start_frame;
    else
        cFrame = p.end_frame;
    end
   
    % Read and display image
    warning off
    im = imread([vPath filesep p.files(cFrame).name]);
    
    % Enhance contrast
    if enhance_contrast
        warning off
        im = adapthisteq(im);
        warning on
    end
    
    figure;
    [p.xlim_c,p.ylim_c] = zoom_image(im,'');
    
    % Select coordinates
    [p.xLine,p.yLine]= choosePoints(im,2,1,'Choose tracking line',...
                                     p.xlim_c,p.ylim_c);
    warning on
    
    % Adjust so that coordinates run from left to right
    if p.xLine(1)>p.xLine(2)
        x_tmp = p.xLine;
        y_tmp = p.yLine;
        p.xLine(2) = x_tmp(1);
        p.xLine(1) = x_tmp(2);
        p.yLine(2) = y_tmp(1);
        p.yLine(1) = y_tmp(2);
        clear x_tmp y_tmp
    end
    
    % Calculate box values
    p.xBox = [p.xLine(1); p.xLine(2); p.xLine(2); p.xLine(1); p.xLine(1)];
    p.yBox = [mean(p.yLine)+ceil(num_pix/2); mean(p.yLine)+ceil(num_pix/2); ...
        mean(p.yLine)-floor(num_pix/2); mean(p.yLine)-floor(num_pix/2); ...
        mean(p.yLine)+ceil(num_pix/2)];
    
    close
    
    save([vPath filesep 'seq_data'],'p')
    
else
    load([vPath filesep 'seq_data.mat'])
    disp(' '); disp('Loaded sequence parameters (seq_data.mat)')
end


%% Interactively adjust smoothing and select peak

if strcmp(p.exp_type,'short')
    cFrame = p.start_frame;
else
    cFrame = p.end_frame;
end

% Read and display image
im = imread([vPath filesep p.files(cFrame).name]);

% Enhance contrast
if enhance_contrast
    warning off
    im = adapthisteq(im);
    warning on
end

% Find raw pixel values
[vals_x,vals] = improfile_horz(im,p.xLine,p.yLine,num_pix);

% Normalize values using preview data
%vals_norm = (vals-min(vals_p(:)))/range(vals_p(:));
vals_norm = (vals-min(vals(:)))/range(vals(:));

% Set initial smoothing value
s.smoothing = smooth_def;

f = figure;

% Selection of smoothing level
while 1
    
    title_txt = 'Adjust smoothing level';
    
    f = display_data(f,im,vals_x,vals_norm,s.smoothing,p,title_txt);
    
    but = questdlg('How is the smoothing level?', ...
                         'Smoothing', 'Less', 'More', 'Fine', 'Fine');
    if strcmp(but,'Less')
        s.smoothing = s.smoothing + 0.1;
    elseif strcmp(but,'More')
        s.smoothing = s.smoothing - 0.1;
    elseif strcmp(but,'Fine')
        break
    else
        return
    end
end

% Choose peak
f = display_data(f,im,vals_x,vals_norm,s.smoothing,p,title_txt);
[xEdge_last,peak] = choose_peak(f,s,p,vals_x,vals_norm);


%% Step thru frames for pipette tip with visual check

f = figure;
set(f,'DoubleBuffer','on')
set(f,'CurrentCharacter','1')

title_txt = {'Left click: opened,   Right click: closed,   p: play', ...
             'space: stop playing,   left arrow: step back'};

% The short movies progress forward, the long ones go backward
if strcmp(p.exp_type,'short')
    frames = p.start_frame:p.end_frame;
else
    frames = p.end_frame:-1:p.start_frame;
end

i = 1;
cFrame = frames(1);

% Loop throgh frames
while 1

    % Read and display image
    im = imread([vPath filesep p.files(cFrame).name]);   
    
    % Enhance contrast
    if enhance_contrast
        warning off
        im = adapthisteq(im);
        warning on
    end

    % Find raw pixel values
    [vals_x,vals] = improfile_horz(im,p.xLine,p.yLine,num_pix);
    
    % Normalize values using preview data
    %vals_norm = (vals-min(vals_p(:)))/range(vals_p(:));
    vals_norm = (vals-min(vals))/range(vals);
    
    [xEdge,peak] = find_peak(vals_x,vals_norm,s.smoothing,xEdge_last,peak);
    
    % Store data from current frame
    s.frame(i) = cFrame;
    s.t(i)     = cFrame./p.frame_rate;
    s.xEdge(i) = xEdge;
    
    % Update plot
    title_txt = {'Left click/right arrow: advance,  p: play', ...
                 'space: stop playing,   left arrow: step back', ...
                 'esc:select manually, return: finish',...
                 ['Frame :' num2str(p.start_frame) '-' num2str(cFrame) '-' ...
                  num2str(p.end_frame)]};
              
    f = display_data(f,im,vals_x,vals_norm,s.smoothing,p,title_txt,xEdge);
   
    
    clear xEdge cFrame
    
    % Play mode: advance frame at regular interval
    if strcmp(get(f,'CurrentCharacter'),'p')
        i   = i + 1;
        pause(.2)
 
    else   
        % Get input from figure window
        [xTmp,yTmp,but] = ginput(1);
        
        % Return
        if isempty(but)
            if i==0
                return   
            else
                break
            end
            
        % Left click/right arrow
        elseif(but==1) || (but==29)
            i  = i + 1;    
            
        % Left arrow key    
        elseif but==28    
            i   = max([1 (i-1)]);
            s.frame = s.frame(1:i);
            s.t     = s.t(1:i);
            s.xEdge = s.xEdge(1:i);
            
        % ESC
        elseif but==27
            disp(''); disp('Select peak manually')
            [s.xEdge(end),peak] = choose_peak(f,s,p,vals_x,vals_norm);
            
        % 'p'
        elseif but==112
            set(f,'CurrentCharacter','p');
            
        % Other
        else
            warning('Not a meanignful key')
        end
        
        clear xTmp yTemp but
    end 
    
    clear im vals title_txt vals_norm vals_x 
    
    % Set for next loop
    xEdge_last  = s.xEdge(end);
    cFrame      = frames(i);
end
    

%% Save data

but = questdlg('Save data?','','Yes','No','Yes');
 
% Reverse order, for long duration recording
if ~strcmp(p.exp_type,'short')
    [tmp,idx] = sort(s.frame);
    
    s.frame = s.frame(idx);
    s.t     = s.t(idx);
    s.xEdge = s.xEdge(idx);
    
    clear idx tmp
end 

% Save
if strcmp(but,'Yes')
    save([vPath filesep 'edge_data'],'s')
else
    return
end

% Clean up
clear xEdge_last
close(f)



function [x_val,peak] = find_peak(vals_x,vals_norm,smoothing,...
                                  x_guess,peak)

% Find splines and derivatives
sp   = spaps(vals_x,vals_norm,10^(-smoothing));
Dsp  = fnder(sp,1);
D2sp = fnder(sp,2);

% Find zeros of 2nd derivative
zero_vals = fnzeros(D2sp,[min(vals_x) max(vals_x)]);

% Cut out second row
zero_vals = zero_vals(1,:);

% Eliminate peaks, if need valleys and visa/versa
if (nargin>4) && peak
    zero_vals = zero_vals(fnval(Dsp,zero_vals)>0);
    
elseif (nargin>4) && ~peak
    zero_vals = zero_vals(fnval(Dsp,zero_vals)<0);
    
end

% % Eliminate zeros in wrong direction
% if (nargin>5) && strcmp(direction,'rightward')
%     zero_vals = zero_vals(zero_vals>x_guess);
%     
% elseif (nargin>5) && strcmp(direction,'leftward')
%     zero_vals = zero_vals(zero_vals<x_guess);
%     
% end

% Difference from guess
diff_val = abs(zero_vals-x_guess);

% Choose zero with smallest difference
x_val = zero_vals(diff_val==min(diff_val));

% Store whether a peak or valley
if fnval(Dsp,x_val)>0
    peak = 1;
else
    peak = 0;
end

function [xlim_c,ylim_c] = zoom_image(im,title_txt)
% Interacively zooms image, returns x and y limits

% Display image
imshow(im)
xlabel('Zoom into your region of interest')
title(title_txt)

% Enter zoom mode
zoom on
pause;

% Determine selected limits
xlim_c = xlim;
ylim_c = ylim;

% Set limits
xlim(xlim_c);
ylim(ylim_c);


function [xLine,vals] = improfile_horz(im,x,y,num_pix)
% Finds mean profile for horizontal measurements by averaging vertically

coord_add = -floor(num_pix/2):ceil(num_pix/2);

for i = 1:length(coord_add)
    [xLine,yLine,val(:,i)] = improfile(im,x,y+coord_add(i));
end

vals = mean(val,2);


function [xEdge,peak] = choose_peak(f,s,p,vals_x,vals_norm)

figure(f)

% Peak/valley selection
while 1
    
    title_txt = '';
 
    subplot(4,1,4)
    title('Select peak or valley to track')
    
    [x_tmp,y_tmp,b_tmp] = ginput(1);
    
    [x_tmp,peak] = find_peak(vals_x,vals_norm,s.smoothing,x_tmp);
    
    subplot(4,1,1:2)
    hold on
    h = plot([x_tmp x_tmp],[min(p.yBox) max(p.yBox)],'g-');
    hold off
    
    subplot(4,1,3)
    hold on
    h(2) = plot([x_tmp x_tmp],ylim,'g-');
    hold off
    
     subplot(4,1,4)
    hold on
    h(3) = plot([x_tmp x_tmp],ylim,'g-');
    hold off
    
    but = questdlg('Is this the correct peak?', ...
                         'Selection','Yes', 'No', 'Yes');
    if strcmp(but,'Yes')
        xEdge = x_tmp;
        break
        
    elseif strcmp(but,'No')
        delete(h)
        
    else
        return
    end
end
 

function f = display_data(f,im,vals_x,vals_norm,smoothing,p,title_txt,xEdge)
    
figure(f);

sp = spaps(vals_x,vals_norm,10^(-smoothing));

% Display frame
subplot(4,1,1:2)
warning off
imshow(im)
title(title_txt)
warning on
xlim(p.xlim_c);
ylim(p.ylim_c);

% Overlay data
hold on
plot(p.xBox,p.yBox,'w-')

if nargin>7
    plot([xEdge xEdge],[min(p.yBox) max(p.yBox)],'g-')
end

hold off

% Plot raw pixel values
subplot(4,1,3)
%     plot(vals_x,vals_norm,'k',...
%          vals_x(idx),valley_func(beta,vals_x(idx)),'r-',...
%          beta(2),beta(3),'or')
plot(vals_x,vals_norm,'k')
hold on
plot(vals_x,fnval(sp,vals_x),'r-')

if nargin>7
    plot([xEdge xEdge],ylim,'g-')
end

hold off
ylabel('Normalized pixel intensity')


subplot(4,1,4)
plot(vals_x,fnval(fnder(sp),vals_x),'r-')
hold on
if nargin>7
    plot([xEdge xEdge],ylim,'g-')
end
hold off
xlabel('Position along line (pix)')
ylabel('First derivative of intensity')
% Plot interpreted values
    



