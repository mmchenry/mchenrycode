function c = inputBehavior(vid_path,c,fishNums)


disp(' ');
disp('Use a video player to target individuals that should be analyzed')
disp('and determine if they are swimming.')
disp(' ');

if isempty(c)
    overWrite = 1;    
else
    overWrite = 0;
end

%% Determine filename for first frame of tifs

% Get filenames
a = dir(vid_path);

% Loop through each filename unitl reaching a tiff image
k = 1;
for i = 1:length(a)
    cName = a(i).name;

    % Define file path from first tif filename with a zero
    % start at 50th frame
    if length(cName)>3 && ...
       strcmp(cName(end-2:end),'tif') && ...
       ~isempty(find(cName=='0',1))
   
        k = k+1;
        if k > 51
            frame_path = [vid_path filesep cName];
            break
        end
        
    elseif i==length(a)
        error('No tif file in given diretcory')
    end
 
end

clear cName k



%% Prompt to choose roi or load from faststart_data

if overWrite
    
    disp(' ')
    disp('Choose region of interest around eyes & trunk for each larva');
    disp(' ========================================================== ');

    fish = 1;
    useCurrent = 1;
    f = figure;
    rois = [];
    im = imread(frame_path);
    imshow(im)
    hold on

    while 1==1
        % Select points from image  
        txt = ['Fish ' num2str(fish) ': Choose region of interest'];
        [roi_x,roi_y] = choosePoints(im,2,2,txt,[],[],useCurrent); 
        
        % Check if no points
        if isempty(roi_x)
            break
        end
        
        % plot rois
        rois = [rois; roi_x(1) roi_x(2) roi_y(1) roi_y(2)];
        hold on
        for i = 1:size(rois,1)
            plot([rois(i,1) rois(i,2) rois(i,2) rois(i,1) rois(i,1)],...
                [rois(i,3) rois(i,3) rois(i,4) rois(i,4) rois(i,3)],'w-')
            h = text(min(rois(i,1:2))+3,min(rois(i,3:4))+7,num2str(i));
            set(h,'Color','w')
        end

        % Store data in c
        c(fish).roi.x = roi_x;
        c(fish).roi.y = roi_y;
        
        % Advance index
        fish = fish + 1;

    end

    gf = getframe(f);
    imwrite(gf.cdata,[vid_path filesep 'regions_of_interest.tif'],'tif')
    
    save([vid_path filesep 'faststart_data'],'c')
    
    clear f h roi_x roi_y useCurrent gf  fish txt rois im 
    
else
    load([vid_path filesep 'faststart_data'])
    
end


%% Get data on individual fish


% Set general parameters
invertImage = 0;
if overWrite
    fishNums  = 1:length(c);
end

% find starting threshold value
roi_x = c(min(fishNums)).roi.x;
roi_y = c(min(fishNums)).roi.y;
im = imread(frame_path);
im_roi = im(round(min(roi_y):max(roi_y)),round(min(roi_x):max(roi_x)));
tVal = graythresh(im_roi);

clear im_roi

for i = fishNums
    
    %Choose threshold level
    roi_x = c(i).roi.x;
    roi_y = c(i).roi.y;
    
    figure;
    set(gcf,'DoubleBuffer','on')
    txt = ['Fish ' num2str(i) ': Choose threshold'];
    [tVal,im_bw] = chooseThresh(im,tVal,invertImage,roi_x,roi_y,txt);
    
    
    % Store tVal in c
    c(i).treshold_level = tVal;
    
    
    % Create binary image (bw) of eyes and trunk
    bw = im2bw(im,tVal);
    if ~invertImage 
        bw = ~bw;
    end
    
    
    % Acquire left eye, find coorindates
    txt = ['Fish ' num2str(i) ': Choose left eye'];
    [left_x,left_y] = choosePoints(im_bw,1,0,txt,roi_x,roi_y);
    bw_tmp  = bwselect(bw,round(left_x),round(left_y));
    L_tmp   = bwlabel(bw_tmp);
    val_tmp = regionprops(L_tmp,'Centroid');
    
    % Store left eye coordinates
    if ~isstruct(val_tmp) && isemtpy(val_tmp)
        error('No spot where you clicked')
    elseif length(val_tmp) > 1 
        error('More than one spot selected')
    else
        c(i).left_eye.x     = val_tmp.Centroid(1);
        c(i).left_eye.y     = val_tmp.Centroid(2);
    end
    
    clear txt bw_tmp L_tmp val_tmp left_x left_y
    
    
    % Acquire right eye, find coorindates
    txt = ['Fish ' num2str(i) ': Choose right eye'];
    [right_x,right_y] = choosePoints(im_bw,1,0,txt,roi_x,roi_y);
    bw_tmp  = bwselect(bw,round(right_x),round(right_y));
    L_tmp   = bwlabel(bw_tmp);
    val_tmp = regionprops(L_tmp,'Centroid');
    
    % Store right eye coordinates
    if ~isstruct(val_tmp) && isemtpy(val_tmp)
        error('No spot where you clicked')
    elseif length(val_tmp) > 1 
        error('More than one spot selected')
    else
        c(i).right_eye.x     = val_tmp.Centroid(1);
        c(i).right_eye.y     = val_tmp.Centroid(2);
    end

    clear txt bw_tmp L_tmp val_tmp right_x right_y
    

    % Acquire trunk, find coorindates
    txt = ['Fish ' num2str(i) ': Choose trunk'];
    [trnk_x,trnk_y] = choosePoints(im_bw,1,0,txt,roi_x,roi_y);
    bw_tmp  = bwselect(bw,round(trnk_x),round(trnk_y));
    L_tmp   = bwlabel(bw_tmp);
    val_tmp = regionprops(L_tmp,'Centroid');
    
    %  Store trunk coordinates
    if ~isstruct(val_tmp) && isemtpy(val_tmp)
        error('No spot where you clicked')
    elseif length(val_tmp) > 1 
        error('More than one spot selected')
    else
        c(i).trnk.x     = val_tmp.Centroid(1);
        c(i).trnk.y     = val_tmp.Centroid(2);
    end
    
    clear txt bw_tmp L_tmp val_tmp trnk_x trnk_y
    
         
    % Prompt user data on this fish
    while 1==1     
        name = ['Fish ' num2str(i) ' data'];
        prompt = {'Was this fish swimming (1 = yes, 0 = no)',...
            'At what frame number did they fast start (0, if never)?'
            };

        answer = inputdlg(prompt,name,1);

        if isempty(answer)
            return
        elseif (~(str2num(answer{1})==0) && ~(str2num(answer{1})==1))
            warning('Answer must be either 0 or 1')
        else
            c(i).swimming               = str2num(answer{1});
            c(i).fast_start_frame_num   = str2num(answer{2});
            break
        end
    end
    
    save([vid_path filesep 'faststart_data'],'c')
    
end





%% function: choosePoints

function [x,y] = choosePoints(im,numLimit,link,title_text,x_roi,y_roi,useCurrent)
%Used for finding coordinate points on a static image 'im'.
% numLimit is the total number of points desired
% link tells how to diplay the points
% link = 0 - unconnected points
% link = 1 - line drawn between 2 most recent points
% link = 2 - bounding rectangle drawn btwn 2 points
% x_roi & y_roi define a region of interest with 2 points
% useCurrent - (1 or 0) to use the current figure window

% Assign defaults for the inputs
if nargin < 7
    useCurrent = 0;
    if nargin < 5
        x_roi = [];
        y_roi = [];
        if nargin < 4
            title_text = ' ';
            if nargin < 3
                link = 0;
                if nargin < 2
                    numLimit = 1000;
                    if nargin < 1
                        error('You need to provide an image');
                    end
                end
            end
        end
    end
end

% Check inputs
if link > 2
    link = 0;
    warning('link cannot be greater than 2, set to 0');
end

if (link==2) && ~(numLimit==2)
    numLimit = 2;
    warning('bounding box requires just 2 coordinates, numLimit set to 2');
end

if nargin == 5 
    error('You need to define both x and y coodinates for your roi');
end

if ~isempty(x_roi) && ...
    (~(length(x_roi)==2) || ~(length(y_roi)==2))
    error('Two coordinates are needed to specify a roi');
end
    

% Display image
if ~useCurrent
    figure;
    if isstruct(im)
        imshow(im.cdata,im.colormap);
    else
        imshow(im)
    end
    if ~isempty(x_roi)
        set(gca,'YLim',[min(y_roi) max(y_roi)])
        set(gca,'XLim',[min(x_roi) max(x_roi)])
    end
    title(title_text)
    set(gcf,'DoubleBuffer','on');
end

% Give instructions
disp(' '); disp(' ');
disp('Left mouse button picks points.');disp(' ');
disp('Right mouse button removes last point.');disp(' ');
disp('Press return when done collecting.')
disp('Press esc to exit');
disp(' '); disp(' ');

% Initiate parameter values
n   = 0;
but = 1;
x = []; y = [];

% Loop through for interactive input
while 1 == 1
    [xi,yi,but] = ginput(1);
    
    if isempty(but)
        break
        
    elseif but==1 % Left click
        
        % Make sure not to exceed limit
        if n+1 > numLimit
            n = numLimit;
        else
            n = n+1;
        end
        
        % Store coordinates
        x(n) = xi;
        y(n) = yi;
        
    elseif but==3 % Right click
        
        % Don't allow n < 0
        if n-1 < 1
            n = 0;
            x = [];
            y = [];
        else
            n = n-1;
            x = x(1:n);
            y = y(1:n);
        end
        
    elseif but==27 % If escape
        x = [];
        y = [];
        break
    end
    
    % Draw data on image
    if isstruct(im)
        imshow(im.cdata,im.colormap);
    else
        imshow(im)
    end
    if ~isempty(x_roi)
        set(gca,'YLim',[min(y_roi) max(y_roi)])
        set(gca,'XLim',[min(x_roi) max(x_roi)])
    end
    title(title_text)
    hold on
    
    if link == 0
        plot(x,y,'+r')
        
    elseif link == 1
        plot(x,y,'o-r')
        
    elseif link == 2
        if length(x)<2
            plot(x,y,'+r')
        else
            plot([x(1) x(2) x(2) x(1) x(1)],...
                 [y(1) y(1) y(2) y(2) y(1)],'r-')
        end
    end
    hold off
    
end
x = x'; y = y';
if ~useCurrent
    close
end

%% function: chooseThresh

function [tVal,img] = chooseThresh(im,tVal,invertImage,x_roi,y_roi,title_txt)
%Used for choosing a threshold value in 'im'.

% Parameters
fillColor  = [.43 .49 1];

% Instructions
disp(' ')
disp('Choose a threshold that fills in the shape of interest')
disp('  as a WHITE blob in a BLACK field, but do not allow');
disp('  the blob to fuse with any other blobs on the image')
disp(' ')
disp('  Up arrow -- increases the threshold')
disp('  Down arrow -- decreases the threshold')
disp('  = key -- micro increase in threshold')
disp('  - key -- micro decrease in threshold')
disp(' ')
disp('  Press return when done')
disp(' ')
    
% Loops each time a button is pushed, until return is pushed
but = 1;

while 1 == 1  
    % convert image to binary at threshold, invert image
    bw = im2bw(im,tVal);
    if ~invertImage
        bw      = ~bw;
    end
   
    % display image
    im_bw        = im(:,:,1);
    im_bw(find(bw)) = 244.*ones(length((find(bw))),1);
    
    cMap = [[0:1/255:1]' [0:1/255:1]' [0:1/255:1]'];
    cMap(245,:) = fillColor;
    imshow(im_bw,cMap);
    if nargin > 3
        set(gca,'YLim',[min(y_roi) max(y_roi)])
        set(gca,'XLim',[min(x_roi) max(x_roi)])
    end
    hold on

    
    
    %title('Make white blob on black field')
    xlabel(['tValue = ' num2str(tVal)])
    title(title_txt)

    %prompt for input, adjust tVal accordingly
    [xi,yi,but] = ginput(1);
    if isempty(but)
        break
    elseif but == 30 %up arrow
        if (tVal + .01) < 1
            tVal = tVal + .01;
        end
    elseif but == 31 %down arrow
        if (tVal - .01) > 0
            tVal = tVal - .01;
        end
    elseif but == 61 %plus sign
        if (tVal + .01) < 1
            tVal = tVal + .001;
        end
    elseif but == 45 %minus sign
        if (tVal - .001) > 0
            tVal = tVal - .001;  
        end
    end
end 

img.cdata     = im_bw;
img.colormap = cMap;

close


%% function: imDisplay


function imBW = imDisplay(im,tVal)
%shows a bitmap image with a threshold level of pixels highlighted in color
imshow(im.cdata,im.colormap);
hold on;

imBW                = ~im2bw(im.cdata,im.colormap,tVal);
im.cdata(find(imBW))= 244.*ones(length((find(imBW))),1);
%cMap                = [[0:1./255:1]' [0:1./255:1]' [0:1./255:1]'];
cmap                = im.colormap;
cmap(245,:)         = fillColor;
imshow(im.cdata,cmap)
hold o
