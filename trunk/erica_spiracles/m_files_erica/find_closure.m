function find_closure(vPath)


%% Parameters

% Number of digits for frame number in filename
num_dig = 4;

% Threshold value (between 0 and 1) for closure cut-off
val_thresh = 0.5;

% Visualize each frame as it is analyzed
vis_frames = true;

% Enahcne image contrast
enhance_con = 1;


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


%% Get info about sequence

answer = inputdlg({'Frame number when opened:',...
                   'Frame number when closed:',...
                   'Threshold value:',...
                   'Frame rate (fps):',...
                   'Correct for jitter (1 for yes, 0 for no)'},'',1,...
                   {'0','16',num2str(val_thresh),'15','1'});

if isempty(answer)
    return
end
               
p.ref.fr_open  = str2num(answer{1}) + 1;
p.ref.fr_close = str2num(answer{2}) + 1;
p.thresh_val   = str2num(answer{3});
p.frame_rate   = str2num(answer{4});
p.jittering    = str2num(answer{5});

% Check input
if isempty(p.ref.fr_open) || isempty(p.ref.fr_close)
    error('You need to provide at least one closed and one opened frame')
end

% Provide instructions
disp(' ')
disp('  First select a point on the cuticle')
disp('  Then a point on in the opening')
disp(' ')
disp('     Right click to erase, return when done')
disp(' ')


%% Select coordinates for pixel profile

% Select line for profile (open frame)
    im = imread_mod([vPath filesep ...
        p.files(p.ref.fr_open==f_nums).name],enhance_con);
    t_txt = ['Opened: ' p.files(p.ref.fr_open==f_nums).name];

    % Prompt for coordinates
    warning off
    figure;
    [xTmp,yTmp]= choosePoints(im,2,1,t_txt);
    close
    warning on

    if length(xTmp) ~= 2
        error('You need to choose two points')
    end

    % xLine and yLine stores the  coordinates 
    p.ref.xLine_open = xTmp';
    p.ref.yLine_open = yTmp';

    clear xTmp yTmp t_txt im

    
% Select line for profile (closed frame)
if ~p.jittering 
    % Just store close coordates as same as open
    p.ref.xLine_close = p.ref.xLine_open;
    p.ref.yLine_close = p.ref.yLine_open;   
    
else
    im = imread_mod([vPath filesep ...
        p.files(p.ref.fr_close==f_nums).name],enhance_con);
    t_txt = ['Closed: ' p.files(p.ref.fr_close==f_nums).name];

    % Prompt for coordinates
    warning off
    figure;
    [xTmp,yTmp]= choosePoints(im,2,1,t_txt);
    close
    warning on

    if length(xTmp) ~= 2
        error('You need to choose two points')
    end
    
    % Find center offset, relative to open frame
    xMean = mean(xTmp) - mean(p.ref.xLine_open);
    yMean = mean(yTmp) - mean(p.ref.yLine_open);
    
    % xLine and yLine values from center offset
    p.ref.xLine_close = p.ref.xLine_open + xMean;
    p.ref.yLine_close = p.ref.yLine_open + yMean;

    clear xTmp yTmp t_txt im xMean yMean
end 

% Select line for profile (first frame)
if ~p.jittering 
    % Just store close coordates as same as open
    p.ref.xLine_first = p.ref.xLine_open;
    p.ref.yLine_first = p.ref.yLine_open;   
    
else
    im = imread_mod([vPath filesep p.files(1).name],enhance_con);
    t_txt = ['First: ' p.files(1).name];

    % Prompt for coordinates
    warning off
    figure;
    [xTmp,yTmp]= choosePoints(im,2,1,t_txt);
    close
    warning on

    if length(xTmp) ~= 2
        error('You need to choose two points')
    end
    
    % Find center offset, relative to open frame
    xMean = mean(xTmp) - mean(p.ref.xLine_open);
    yMean = mean(yTmp) - mean(p.ref.yLine_open);
    
    % xLine and yLine values from center offset
    p.ref.xLine_first = p.ref.xLine_open + xMean;
    p.ref.yLine_first = p.ref.yLine_open + yMean;

    clear xTmp yTmp t_txt im xMean yMean
end 


%% Get profiles of openings and closures

% Closure profile, stored in 'vals'
    im_cl = imread_mod([vPath filesep ...
               p.files(p.ref.fr_close==f_nums).name],enhance_con);
    p.ref.vals.close = find_improfile(im_cl,p.ref.xLine_close,p.ref.yLine_close);

% Opening profile, stored in 'vals'
    im_op = imread_mod([vPath filesep ...
               p.files(p.ref.fr_open==f_nums).name],enhance_con);
    p.ref.vals.open = find_improfile(im_op,p.ref.xLine_open,p.ref.yLine_open);
    
% Store range of value difference btwn opened and closed
    p.ref.vals.range = max(abs(p.ref.vals.close - p.ref.vals.open));
    
    clear im_cl im_op



%% Step through frames for acquisition

f = figure;
set(f,'DoubleBuffer','on')
set(f,'CurrentCharacter','1')

% Grab 'closed' frame
im_close = imread_mod([vPath filesep ...
                p.files(p.ref.fr_close==f_nums).name],enhance_con);

% Grab initial frame
%im0 = imread([vPath filesep p.files(1).name]);
%first_frame = find(p.ref.fr_open(1)==f_nums,1,'first');
%last_frame  = 

%im0 = imread([vPath filesep p.files(1).name]);

% Set 'last' image
im_last = imread_mod([vPath filesep p.files(1).name],enhance_con);


% Initialize variables
xOff = 0;
yOff = 0;
Disps = [];

cFrame = 1;
data_idx = 1;
    
while 1
    
    % Read current frame
    im = imread_mod([vPath filesep p.files(cFrame).name],enhance_con);

    % Correct for image movement
    if p.jittering
        
        % Find displacement
        [xDisp,yDisp] = imshift(im,im_last); 
        
        % Calcuate resultant
        cDisp = sqrt(xDisp^2 + yDisp^2);
        
        % Check that displacement isn't out of bounds
        if ((1-data_idx) > 50) && (cDisp > 1.5*max(Disp))
            xDisp = 0;
            yDisp = 0;
            disp('Skipping on displacement . . .')
            
        % Otherwise, store away typical max displacements for first 50 frames
        else
            Disps = [Disps; cDisp];
            
        end
                
        xOff = xOff + xDisp;
        yOff = yOff + yDisp;
        
    else
        xOff = 0;
        yOff = 0;
    end
    
    % Find pixel values
    val = find_improfile(im,p.ref.xLine_first,p.ref.yLine_first);
    
    % get frame number
    fr_num = str2num(p.files(cFrame).name(end-3-num_dig:end-3));
    
    % Normalize pixel values 
    val_norm = abs(val - p.ref.vals.open)./p.ref.vals.range;

    % Binary on whether closed
    closed = mean(val_norm) > p.thresh_val;
    
    
    % Visualize current data ----------------------------------------------
    subplot(3,3,1:6)
    warning off
    imshow(im)
    
    warning on
    hold on
    h1 = plot(p.ref.xLine_first+xOff,...
        p.ref.yLine_first+yOff,'r-');
    hold off
    title([p.files(cFrame).name '       Click any key to pause'])
    set(h1,'LineWidth',3)
    
    
    subplot(3,3,7:8)
    h3 = plot(val_norm,'r');
    ylabel('Pixel value')
    xlabel('Position (pixels)')
    
    %         title(['Mean val_norm = ' ...
    %                  num2str(mean(p.profile(i).val_norm)) ])
    set(h3,'LineWidth',2)
    ylim([-.1 1.2])
    
    subplot(3,3,9)
    h4 = bar(mean(val_norm),'r');
    title('Mean pixel value')
    ylim([-.1 1.2])
    hold on
    plot([.5 1.5],[p.thresh_val p.thresh_val],'k-')
    hold off
    
    if closed
        subplot(3,3,1:6)
        hold on
        h2 = plot([1 size(im,2)-1 size(im,2)-1 1 1],...
            [1 1 size(im,1) size(im,1) 1],'g-');
        hold off
        set(h2,'LineWidth',8)
        set(h1,'Color','g')
        set(h3,'Color','g')
        set(h4,'FaceColor','g')
        
    end    
    % End of visualization ------------------------------------------------
    
    
    % Pause to allow key press
    if i==1
        pause(.5)
    else
        pause(.2)
    end
    
    % Break loop, if key is pressed
    if ~strcmp(get(f,'CurrentCharacter'),'1')
        
        % Recreate figure window
        close
        f = figure;
        set(f,'DoubleBuffer','on')
        
        % number of frames to jump back
        numJump = 3;
        
        % Jump back a few frames
        cFrame    = cFrame - numJump;
        data_idx  = data_idx - numJump; 
        beep
        disp(['Jumping back by ' num2str(numJump) ' frames']);
        
        % Read frame
        im = imread_mod([vPath filesep p.files(cFrame).name],enhance_con);
        t_txt = ['Reset position: ' p.files(cFrame).name];
        
        % Prompt for coordinates
        warning off
        subplot(3,3,1:6)
        [xTmp,yTmp]= choosePoints(im,2,1,t_txt);
        warning on
        
        if length(xTmp) ~= 2
            error('You need to choose two points')
        end
        
        % Find center offset, relative to open frame
        xMean = mean(xTmp) - mean(p.ref.xLine_open);
        yMean = mean(yTmp) - mean(p.ref.yLine_open);
        
        % xLine and yLine values from center offset
        p.ref.xLine_first = p.ref.xLine_open + xMean;
        p.ref.yLine_first = p.ref.yLine_open + yMean;
        
        % Reset offset to zero
        xOff = 0;
        yOff = 0;
        
        % Prompt for new threshold value
        answer = inputdlg('Threshold value:','',1,{num2str(p.thresh_val)});
        if isempty(answer)
            return
        else
           p.thresh_val =  str2num(answer{1});
        end
        
        % Prep for next loop
        im_last = im;
        set(f,'CurrentCharacter','1')

        clear xTmp yTmp t_txt im xMean yMean numJump answer
    
    % Store values, if no key pressed
    else
        
        % Store raw pixel values
        p.profile(data_idx).vals = val;
        
        % Normalize pixel values (for plotting)
        p.profile(data_idx).val_norm = val_norm;
        
        % Store binary on whether closed
        p.time(data_idx)   = fr_num/p.frame_rate;
        p.closed(data_idx) = closed;
        
        
        % Prep for next loop
        im_last = im;
        cFrame  = cFrame + 1;
        data_idx = data_idx + 1;
        
    end
   
    clear im val val_norm closed
    
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

function [x,y] = imshift(im,im0)
% Measures the displacement of image 'im', relative to
% reference image 'im0'

% Parameters
maxArea  = 50;

% % Enhance contrast of images
% warning off
% im   = adapthisteq(im);
% im0  = adapthisteq(im0);
% warning on

% Eliminate small particles
warning off
se = strel('disk',5);
im_tmp = imopen(im,se);
im0_tmp = imopen(im0,se);
warning on

% Find threshold level
level = max([.3 1.5*graythresh(im0)]);

% Threshold
bw = im2bw(im_tmp,level);
bw0 = im2bw(im0_tmp,level);

% Gather properties
warning off
stats = regionprops(bw,'Area','Eccentricity','Centroid');
stats0 = regionprops(bw0,'Area','Eccentricity','Centroid');
warning on

% Make properties into vectors
for i = 1:length(stats)
    areas(i,1) = stats(i).Area;
    eccs(i,1)  = stats(i).Eccentricity;
    Xs(i,1)    = stats(i).Centroid(1);
    Ys(i,1)    = stats(i).Centroid(2);
end

for i = 1:length(stats0)
    areas0(i,1) = stats0(i).Area;
    eccs0(i,1)  = stats0(i).Eccentricity;
    Xs0(i,1)    = stats0(i).Centroid(1);
    Ys0(i,1)    = stats0(i).Centroid(2);
end

% Check for enough blobs
if (length(areas) < 3) || (length(areas0) < 3) 
    error({'not enough features in image to correct for jitter,',...
            'try a lower threshold'})
end

% Meshgrid the properties (normalized)
[A1,A2] = meshgrid(areas./max(areas),areas0./max(areas0));
[E1,E2] = meshgrid(eccs./max(eccs),eccs0./max(eccs0));

% Define index of differences, where the rows correspond to 
% elements in im0
diff_idx = abs(A1-A2).*abs(E1-E2);

% Step through elements in im0
for i = 1:size(diff_idx,1)
    
    % The element in im that matches current element in im0
    best_col = find(diff_idx(i,:)==min(diff_idx(i,:)),1,'first');
    
    x_diff(i,1) = Xs(best_col) - Xs0(i);
    y_diff(i,1) = Ys(best_col) - Ys0(i);
end

% Find lowest half of values
disp = sqrt(x_diff.^2 + y_diff.^2);
[disp,idx] = sort(disp);
idx = idx(1:round(length(disp)/4));

% Take displacement as mean of lowest half of values
x = mean(x_diff(idx));
y = mean(y_diff(idx));



function vals = find_improfile(im,x,y)
% Finds mean profile around specified coordinates

coord_add = -1:1;

for i = 1:length(coord_add)
    val1(:,i) = improfile(im,x+coord_add(i),y+coord_add(i));
end

for i = 1:length(coord_add)
    val2(:,i) = improfile(im,x-coord_add(i),y-coord_add(i));
end

vals = mean([val1 val2],2);