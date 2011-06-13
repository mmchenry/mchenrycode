function choose_closure(vPath,frame_rate)
% Used in place of find_closure, in case of challenging movie.



%% Parameters

% Number of digits for frame number in filename
num_dig = 4;

% Threshold value (between 0 and 1) for closure cut-off
val_thresh = 0.5;

% Visualize each frame as it is analyzed
vis_frames = true;

% Enhance image contrast
enhance_con = 1;

% Prompt for frame rate
if nargin<2
    answer = inputdlg({'Frame rate (fps):'},'',1,{'15'});
    
    if isempty(answer)
        return
    end
    
    p.frame_rate   = str2num(answer{1});
    
else
    p.frame_rate = frame_rate;
end


%% Select path

% Prompt for directory
if nargin < 1
   [vPath] = uigetdir(pwd,'Choose video directory');
end

% Prompt if data file is present
if ~isempty(dir([vPath filesep 'profile_data.mat']))
    button = questdlg('Data file already exists.  Overwrite?','Warning!',...
                      'Yes','Cancel','Cancel');
    if strcmp(button,'Yes')
        delete([vPath filesep 'profile_data.mat'])
    else
        return
    end
    
    clear button
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


%% Make vector of frame numbers

for i = 1:length(p.files)
    f_nums(i) = 1+str2num(p.files(i).name(end-4-num_dig+1:end-4));
end


%% Step through frames for acquisition

f = figure;
set(f,'DoubleBuffer','on')
set(f,'CurrentCharacter','1')

title_txt = {'Left click: opened,   Right click: closed,   p: play', ...
             'space: stop playing,   left arrow: step back'};

play_vid = 0;

% Initialize variables
xOff = 0;
yOff = 0;
Disps = [];

cFrame = 1;
closed = 0;
p.closed = [];
    
while 1
    
    % Read current frame
    im = imread_mod([vPath filesep p.files(cFrame).name],enhance_con);
    
    % get frame number
    fr_num = str2num(p.files(cFrame).name(end-3-num_dig:end-3));    
    
    % Show video frame 
    warning off
    imshow(im)
    title(title_txt)
    xlabel([p.files(cFrame).name])
    warning on
    
    % Draw boarder 
    hold on
    h1 = plot([1 size(im,2)-1 size(im,2)-1 1 1],...
            [1 1 size(im,1) size(im,1) 1],'w-');
    hold off
    set(h1,'LineWidth',8)
    
    % Draw boarder, if backing thru frames
    if (cFrame < length(p.closed)) && p.closed(cFrame)
        set(h1,'Color','g')
        
    elseif (cFrame < length(p.closed)) && ~p.closed(cFrame)
        set(h1,'Color','r')
        
    end
    
    %set(f,'Selected','on')
    figure(f)    
    
    if strcmp(get(f,'CurrentCharacter'),'p')
        closed   = 0;
        cFrame   = cFrame + 1;

        set(h1,'Color','r')
        
    else   
        % Get input from figure window
        [xTmp,yTmp,but] = ginput(1);
        
        % Left click
        if but==1
            closed = 0;
            cFrame  = cFrame + 1;
            
        % Right click
        elseif but==3
            closed = 1;
            cFrame  = cFrame + 1;
            
        % Left arrow key    
        elseif but==28    
            cFrame   = cFrame - 1;
            closed   = p.closed(cFrame);
            
        % ESC
        elseif but==27
            return
            
        % 'p'
        elseif but==112
            closed = 0;
            set(f,'CurrentCharacter','p');
            
        % Other
        else
            closed = 0;
            warning('Not a meanignful key')
        end
        
        clear xTmp yTemp but
    end
    
    % Update frame color
    if closed
        set(h1,'Color','g')
    else
        set(h1,'Color','r')
    end
        
    pause(.2)
   
    
    % Store binary on whether closed
    p.time(cFrame)   = fr_num/p.frame_rate;
    p.closed(cFrame) = closed;  
    

   % set(f,'CurrentCharacter','1')

    clear im val val_norm closed
        
    if cFrame > max(f_nums)
        break
    end
end


%% Save results

save([vPath filesep 'profile_data'],'p')

disp('Done !')
close


%% Visualize

see_closure(vPath)

return

function im = imread_mod(im_path,enhance_con)
% Modified version of im_read

% Read image
im = imread(im_path);

% Enhance contrast of image
if enhance_con
    warning off
    im   = adapthisteq(im);
    warning on
end
