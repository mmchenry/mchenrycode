function find_closure(vPath)


%% Parameters

% Number of digits for frame number in filename
num_dig = 4;

% Threshold value (between 0 and 1) for closure cut-off
val_thresh = 0.5;

% Visualize each frame as it is analyzed
vis_frames = true;


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


%% Get coordinates from opening frames

answer = inputdlg({'Frame number when opened:',...
                   'Frame number when closed:',...
                   'Threshold value:',...
                   'Frame rate (fps):'},'',1,...
                   {'75','79',num2str(val_thresh),'15'});

p.ref.fr_open  = str2num(answer{1});
p.ref.fr_close = str2num(answer{2});
p.thresh_val   = str2num(answer{3});
p.frame_rate   = str2num(answer{4});

% Check input
if isempty(p.ref.fr_open) || isempty(p.ref.fr_close)
    error('You need to provide at least one closed and one opened frame')
    
elseif length(p.ref.fr_open) ~= length(p.ref.fr_close)
    error('You need to provide an equal number of closed and opened frames')
    
end

% Provide instructions
disp(' ')
disp('  First select a point on the cuticle')
disp('  Then a point on in the opening')
disp(' ')
disp('     Right click to erase, return when done')
disp(' ')

% Loop though opening frames to select line for profile
for i = 1:length(p.ref.fr_open)
    im = imread([vPath filesep ...
        p.files(p.ref.fr_open(i)==f_nums).name]);
    t_txt = [p.files(p.ref.fr_open(i)==f_nums).name];
    
    % Prompt for coordinates
    warning off
    [xTmp,yTmp]= choosePoints(im,2,1,t_txt);
    warning on
    
    if length(xTmp) ~= 2
        error('You need to choose two points')
    end
    
    % xLine and yLine store the two coordinates for each 'open' frame
    p.ref.xLine(i,:) = xTmp';
    p.ref.yLine(i,:) = yTmp';
    
    clear xTmp yTmp t_txt im
end





%% Get profiles of openings and closures

% Closure and opening profiles, stored in 'vals'
for i = 1:length(p.ref.fr_close)
    im_cl = imread([vPath filesep ...
        p.files(p.ref.fr_close(i)==f_nums).name]);
    p.ref.vals(i).close = find_improfile(im_cl,p.ref.xLine(i,:),p.ref.yLine(i,:));
    
    im_op = imread([vPath filesep ...
        p.files(p.ref.fr_open(i)==f_nums).name]);
    p.ref.vals(i).open = find_improfile(im_op,p.ref.xLine(i,:),p.ref.yLine(i,:));
    
    % Store range of value difference btwn opened and closed
    p.ref.vals(i).range = max(abs(p.ref.vals(i).close - p.ref.vals(i).open));
    
    clear im_cl im_op
end


%% Step through to capture images

j = 1;

if vis_frames
    figure
    set(gcf,'DoubleBuffer','on')
    pause(.1)
end

% Grab 'closed' frame
im_close = imread([vPath filesep ...
    p.files(p.ref.fr_close(j)==f_nums).name]);

for i = 1:length(p.files)
    
    % Read current frame
    im = imread([vPath filesep p.files(i).name]);
    
    % Find pixel values
    val = find_improfile(im,p.ref.xLine(j,:),p.ref.yLine(j,:));
    
    % get frame number
    fr_num = str2num(p.files(i).name(end-3-num_dig:end-3));
    
    % Store raw pixel values
    p.profile(i).vals = val;
    
    % Normalize pixel values (for plotting)
    p.profile(i).val_norm = abs(val- p.ref.vals(j).open)./...
                            p.ref.vals(j).range;

    % Store binary on whether closed
    p.time(i)   = fr_num/p.frame_rate;
    p.closed(i) = mean(p.profile(i).val_norm) > p.thresh_val;
    
    % Visualize analysis
    if vis_frames
        
        subplot(3,3,1:6)
        warning off
        imshow(im)
        
        warning on
        hold on 
        h1 = plot(p.ref.xLine(j,:),p.ref.yLine(j,:),'r-');
        hold off
        title([p.files(i).name])
        set(h1,'LineWidth',3)

        
        subplot(3,3,7:8)
        h3 = plot(p.profile(i).val_norm,'r');
        ylabel('Pixel value')
        xlabel('Position (pixels)')
        
%         title(['Mean val_norm = ' ...
%                  num2str(mean(p.profile(i).val_norm)) ])
        set(h3,'LineWidth',2)
        ylim([-.1 1.2])
        
        subplot(3,3,9)
        h4 = bar(mean(p.profile(i).val_norm),'r');
        title('Mean pixel value')
        ylim([-.1 1.2])
        hold on
        plot([.5 1.5],[p.thresh_val p.thresh_val],'k-')
        hold off
        
        if p.closed(i)
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
        
        
        if i==1           
            pause(.5)
        else 
            pause(.01)
        end
                
    else
        % Update status
        disp(['Done ' num2str(i) ' of ' num2str(length(p.files))])
    end
    
   
    clear im val
    
end

%% Save results

save([vPath filesep 'profile_data'],'p')

disp('Done !')
close


%% Visualize

see_closure(vPath)

return




function vals = find_improfile(im,x,y)


coord_add = -1:1;

for i = 1:length(coord_add)
    val1(:,i) = improfile(im,x+coord_add(i),y+coord_add(i));
end

for i = 1:length(coord_add)
    val2(:,i) = improfile(im,x-coord_add(i),y-coord_add(i));
end

vals = mean([val1 val2],2);