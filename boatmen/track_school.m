function track_school
% Tracks the position and orientation of a group of fish


%% Parameter values

% Header for image filenames
nameHead = 'CIMG';

% Extension for image files
nameSuffix = 'JPG';

% Visualize acquisition
visSteps = 1;

% Max number of frames for creating the mean image
maxFrames = 1000;


%% Define paths

% Path to this m-file
mPath = '/Volumes/Docs/Projects/rheotaxis/m_files';

% Path to video frames 
vPath = '/Volumes/Docs/Projects/rheotaxis/treated';

% Load filenames for frames
a = dir([vPath filesep  nameHead '*.' nameSuffix]);


%% Create figure window

f = figure;
set(f,'DoubleBuffer','on');


%% Interactively define roi

% Look for roi data
a2 = dir([vPath filesep 'roi.mat']);

if isempty(a2)
    
    % Show first frame in figure window
    im = rgb2gray(imread([vPath filesep a(1).name]));
    warning off
    imshow(im)
    warning on
    hold on
    
    % Select point 1
    title('Choose 2 points to define roi')
    [x,y] = ginput(1);
    h = plot(x,y,'r+');
    
    % Select point 2
    [x(2,1),y(2,1)] = ginput(1);
    
    % Define roi
    xVals = [min(round(x)) max(round(x))];
    yVals = [min(round(y)) max(round(y))];
    
    % Overlay roi
    delete(h)
    h = plot([xVals(1) xVals(2) xVals(2) xVals(1) xVals(1)],...
        [yVals(2) yVals(2) yVals(1) yVals(1) yVals(2)],'r-');
    pause(0.5)
    hold off
    
    % Store data
    roi.xVals = xVals;
    roi.yVals = yVals;
    
    % Save roi data
    save([vPath filesep 'roi.mat'],'roi')
    
    % Clear variables
    clear h x y xVals yVals im a2 roi
    
end


%% Create  or load mean image

% Look for mean image
a2 = dir([vPath filesep 'meanImage.tif']);

% Calculate mean image does not exist
if isempty(a2)   
    
    % Define list of frame numbers, depending on max number of frames
    % requested
    if length(a) > maxFrames
        dframe = round(length(a)/maxFrames);
        frames = 1:dframe:length(a);
        clear dframe
    else
        frames = 1:length(a);
    end
    
    % Create waitbar
    h = waitbar(0,...
            ['Mean image, Frame: ' num2str(frames(1))],...
             'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    % Create sum image based on first frame
    [imCurr,tmp] = imread([vPath filesep a(frames(1)).name]);
    imSum = double(imCurr);
    clear imCurr tmp
    
    % Loop through frames 
    for i = 2:length(frames)
        
        % Add current frame to sum image
        [imCurr,tmp] = imread([vPath filesep a(frames(i)).name]);
        imSum        = imSum + double(imCurr);
        clear tmp imCurr
        
        % Update status bar
        h = waitbar(i/length(frames),h,...
            ['Mean image, Frame: ' num2str(i) ' of ' num2str(length(frames))]);
        
        % Quit m-file, if cancel button pushed
        if getappdata(h,'canceling')
            close force
            return
        end
        
    end
    
    % Calculate mean from sum image
    imMean = uint8(round(imSum./length(frames)));
    
    imMean = imMean(:,:,1);
    
    % Write image to movie dir
    imwrite(imMean,[vPath filesep 'meanImage.tif'],'tif',...
            'Compression','none');
    
    close force
    clear frames h i imSum
        
    %imMean = rgb2gray(imMean);
      
    
% Load mean image, if present
else
    
    disp(' ')
    disp('Loading mean image . . .');
    imMean = imread([vPath filesep 'meanImage.tif']);
    
end


%% Interactively define other areas, threshold

% Load roi
load([vPath filesep 'roi.mat'])

% Look for parameters values, p
a2 = dir([vPath filesep 'params.mat']);

if isempty(a2)

    % Load default variables
    load([mPath filesep 'params.mat'])
    
    % Load first frame of video
    im = rgb2gray(imread([vPath filesep a(1).name]));
    
    % Prompt user to modify parameter values
    an = ' ';
    while ~strcmp(an,'Yes: good image')
        
        % Display image
        visImage(im,imMean,roi,p,1);
        
        % Prompt for input
        an = questdlg('Okay to proceed?','Value change?',...
            'Yes: good image','No: change params','Cancel',...
            'Yes: good image');
        
        % Quit, if cancel
        if strcmp(an,'Cancel')
            
            return
            
        % Prompt for param values, update p
        elseif strcmp(an,'No: change params')
            prompt = {'Min area (pix)','Max area (pix)',...
                      'Threshold value (0 to 1)'};
            pVals = {num2str(p.minArea),num2str(p.maxArea),num2str(p.tVal)};
            an2 = inputdlg(prompt,'Change parameter values',1,pVals);
            
            p.minArea = str2num(an2{1});
            p.maxArea = str2num(an2{2});
            p.tVal = str2num(an2{3});
            
            clear prompt pVals an2
        end
    end 
    
    % Values for sequence
    save([vPath filesep 'params.mat'],'p')
    
    % Default values
    save([mPath filesep 'params.mat'],'p')
    
    % Clear parameter values
    clear an im p

end


%% Step through frames 

% Load roi
load([vPath filesep 'roi.mat'])

% Load parameter values, p
load([vPath filesep 'params.mat'])

% Loop through frames
for i = 1:length(a)
   
    % Load image
    im = rgb2gray(imread([vPath filesep a(i).name]));
   
    % Get image stats
    [im,cmap,bw,stats] = visImage(im,imMean,roi,p,0);
    
    % Step through each blob
    if visSteps
        warning off
        imshow(im);
        hold on
        warning on
    end
    
    % Loop through each blob
    for j = 1:length(stats)
        
        % Center of area of blob
        xC = stats(j).Centroid(1);
        yC = stats(j).Centroid(2);
        
        % Coordinates of weighted (by pixel intensity) centroid
        xW = stats(j).WeightedCentroid(1);
        yW = stats(j).WeightedCentroid(2);
        
        % Define unit vector for orientation
        unitVect = [xW-xC yW-yC]./norm([xC-xW yC-yW]);
        
        % Approximate the nose position from the major axis length
        xTail = unitVect(1) * 0.9*(stats(j).MajorAxisLength/2) + xW;
        yTail = unitVect(2) * 0.9*(stats(j).MajorAxisLength/2) + yW;
        
        % Plot the results
        if visSteps
            h = plot(xTail,yTail,'ro');
            set(h,'MarkerfaceColor','r');
            hold on
            h = plot([xC xTail],[yC yTail],'r-');
            set(h,'LineWidth',2)
            title(['Frame ' num2str(i) ' of ' num2str(length(a))])
        end
    end
    
    % Pause, if visualizing acquisition
    if visSteps
        pause(.2)
        hold off
    else
        disp(['Done ' num2str(i) ' of ' num2str(length(a))])
    end
    
end


%% TO DO:Post-processing




function [im,cmap,bw2,stats] = visImage(im,imMean,roi,p,visOn)


% TO DO: Add image subtraction


% Crop image
im = im([roi.yVals(1):roi.yVals(2)],[roi.xVals(1):roi.xVals(2)],:);

% Inverted image used to find weighted centroid
imInvert = imcomplement(im);

% Crop mean image
imMean = imMean([roi.yVals(1):roi.yVals(2)],[roi.xVals(1):roi.xVals(2)],:);

% Substract mean image
warning off
im     = imcomplement(imsubtract(imMean,im));
warning on

% Binary image for thresholding
bw1 = ~im2bw(im,p.tVal);

% Fill holes
%bw1 = imfill(bw1,'holes');

% Binary image for selecting areas
bw2 = bw1;

% Identify blobs
%L = bwlabel(bw2,8);
stats = regionprops(bw2,'Area','PixelList');

% Area variables
minA = 100000;
maxA = 0;

% Step trhough blobs, reject based on area
for i = 1:length(stats)
    
    if (stats(i).Area < p.minArea) ||  (stats(i).Area > p.maxArea)
        %idx = L(:)==i;
        %bw2(L==i) = 0;
        bw2(stats(i).PixelList(:,2),stats(i).PixelList(:,1)) = 0;
    else
        minA = min([minA stats(i).Area]);
        maxA = max([maxA stats(i).Area]);
    end
        
end

% Redefine stats, after rejecting blobs by area
clear stats
stats = regionprops(bw2,imInvert,'Area','Centroid','Orientation',...
                    'WeightedCentroid','MajorAxisLength');

% Define pixel values, based on threshold and area
im(bw1(:)) = 213;
im(bw2(:)) = 234;

% Alter colormap
cmap        = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
cmap(235,:) = [1 0 0];
cmap(214,:) = [255 180 153]./255;
    
% Enhance contrast of image 
im = imadjust(im);

% Display and redefine colormap
if visOn
    warning off
    imshow(im)
    warning on   
    colormap(cmap)
    title(['Min blob area = ' num2str(minA) '    Max blob area = ' num2str(maxA)])
end



return
% Load subtraction image
if lt_on
    imSub  = imread([mov.dirPath filesep,'meanImage_lt.tif']);   
else
    imSub  = imread([mov.dirPath filesep,'meanImage_dk.tif']);
end

% Adjust grayscale values and convert to double
im     = img.cdata;
warning off
im     = imsubtract(imMean,im);
warning on

if invert
    im = imcomplement(im);
end

if min(im(:))==255
    warning(['Frame ' num2str(fNum) ' is completely white']);
else
    im = imadjust(im,[double(min(im(:)))/255;255/255],[10/255;255/255]);
end
img.cdata = im;


