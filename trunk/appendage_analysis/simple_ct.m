function simple_ct(f_name,imPath)
% CT analysis (that does not include an interactive component) to calculate the moment of inertia


%% Parameters

num_digits = 4;     % Digits at end of image filenames


%% Get info about image files

% Prompt for first slice, if not given
if nargin < 2  
    [f_name,imPath,fIndex] = uigetfile({'*.jpg';'*.tif'},...
        'Choose first image in sequence');
    if ~fIndex
        return
    end
end

% Determine number of slices, file information
[tPath,tName,tExt,tVers] = fileparts([imPath filesep f_name]);
prefix = tName(1:end-num_digits);
dirOutput = dir([imPath filesep prefix '*' tExt]);
fileNames = {dirOutput.name}';
numSlices = numel(fileNames);

clear tPath tName tExt tVers dirOutput prefix num_digits f_name

% Check for images
if isempty(fileNames)
    error(['No image files match ' prefix '*' tExt ' in ' imPath])
end


%% Define stack info

% Look for data file
a = dir([imPath filesep 'stackData.mat']);

if isempty(a)
    
    % Prompt user for input
    prompt = {'Last empty image number',...
              'Spacing between z-slices (mm)',...
              'XY calibration constant (mm/pix)',...
              'Frame num with medial joint',...
              'Frame num with lateral joint',...
              'Apply masks to images (1 or 0)?'};
    defaults = {'1','.0123','.0123','','','0'};
    answer = inputdlg(prompt,'Input stack data',1,defaults);
 
    if isempty(answer{1})
        return
    end
    
    % Store answers, clear unnecessary variables
    s.startSlice    = str2num(answer{1})+1;
    s.endSlice      = numSlices;
    s.zSliceSpacing = str2num(answer{2}).*10^-3;
    s.XYcalconst    = str2num(answer{3}).*10^-3;
    s.frMedJoint    = str2num(answer{4});
    s.frLatJoint    = str2num(answer{5});
    s.useMasks      = str2num(answer{6});   
    
    clear prompt answer defaults a
    
    % Read first image, convert to grayscale
    im = imread([imPath filesep fileNames{s.startSlice}]);
    if size(im,3)==3
        im = rgb2gray(im);
    end
    
    % Define image dimensions
    s.imWidth   = size(im,2);
    s.imHeight  = size(im,1);  
    
    % 24-bit images require a mask and filtering of the gray field
    if s.useMasks
        
        % Check field pixel value
        figure;       
        imDisplay(im)
        hold on
        plot(round(s.imWidth/4),round(s.imHeight/4),'ro')
        title('Position of gray pixel (return to continue)')
        pause
        close
        
        s.grayVal = im(round(s.imWidth/4),round(s.imHeight/4));
        
        % Select masks for static text on slices
        disp(' ')
        disp('Select top left mask (1 of 2) ====================================')
        [s.rMsk1,s.cMsk1] = makeMask_sqr(im,1);
        
        disp(' ')
        disp('Select bottom right mask (2 of 2) ====================================')
        [s.rMsk2,s.cMsk2] = makeMask_sqr(im,3);
        
        % Check that masks were provided
        if isempty(s.rMsk1) || isempty(s.rMsk2)
            error('Mask(s) skipped!')
        end
    end
    
    clear im
    
    % Select point for medial joint
    im = imread([imPath filesep fileNames{s.frMedJoint}]);
    if size(im,3)==3
        im = rgb2gray(im);
    end
    f = figure;
    set(f,'DoubleBuffer','on')
    imDisplay(im)
    title('Zoom, then select white dot for joint position')
    hold on
    zoom on
    [x,y,b] = ginput(1);
    plot(x,y,'r+');
    pause(1)
    close
    
    s.xMedJoint = x;
    s.yMedJoint = y;
    
    clear x y b im
    
    
    % Select point for lateral joint
    im = imread([imPath filesep fileNames{s.frMedJoint}]);
    if size(im,3)==3
        im = rgb2gray(im);
    end
    f = figure;
    set(f,'DoubleBuffer','on')
    imDisplay(im)
    title('Zoom, then select white dot for joint position')
    hold on
    zoom on
    [x,y,b] = ginput(1);
    plot(x,y,'r+');
    pause(1)
    close
    
    s.xLatJoint = x;
    s.yLatJoint = y;
    
    clear x y b im
    
    
    % Save stack data
    disp(' '); disp('Saving stackData.mat')
    save([imPath filesep 'stackData.mat'],'s')
    
    % Clear variables
    clear s
end
    
clear numSlices


%% Acquire pixel values

% Check for existing data file
a = dir([imPath filesep 'I_by_slice.mat']);

if isempty(a)
    
    % Load stack data
    load([imPath filesep 'stackData.mat'])
    
    % Slice list
    sliceList = s.startSlice:s.endSlice;
    
    % Unit volume
    dV = s.XYcalconst^2 .* s.zSliceSpacing;
    
    % Pixel coordinates, offset to voxel center (pix)
    [X,Y] = meshgrid(1:s.imWidth,1:s.imHeight);
    Z     = sliceList - 0.5;
    X     = X - 0.5;
    Y     = Y - 0.5;
    
    % Extract coordinates for joint (pix)
    Jx = mean([s.xLatJoint  s.xMedJoint]);
    Jy = mean([s.yLatJoint  s.yMedJoint]);
    Jz = mean([s.frLatJoint s.frMedJoint]);
    
    % Convert linear dimensions to m
    X  = X  .* s.XYcalconst;
    Y  = Y  .* s.XYcalconst;
    Z  = Z  .* s.zSliceSpacing;
    Jx = Jx .* s.XYcalconst;
    Jy = Jy .* s.XYcalconst;
    Jz = Jz .* s.zSliceSpacing;
    
    % Start timer
    tic
      
    % Loop through slices
    for i = 1:length(sliceList) 
        
        sNum = sliceList(i);
        
        % Read image, convert to grayscale
        im = imread([imPath filesep fileNames{sNum}]);
        if size(im,3)==3
            im = rgb2gray(im);
        end
        
        if s.useMasks
            % Fill masks
            im(s.rMsk1,s.cMsk1) = s.grayVal;
            im(s.rMsk2,s.cMsk2) = s.grayVal;
            
            % Replace light gray for black
            im(im==s.grayVal) = 0;
            
            % Filter noise at edge of morphology
            se = strel('disk',1);
            im = imerode(im,se);
        end 
        
        % Radius (wrt joint) for each pixel in image
        r = ( (X-Jx).^2 + (Y-Jy).^2 + (Z(i)-Jz).^2 ).^0.5;
        
        % Product of radius^2, dV and pixel intensity
        prod1 = r.^2 .* double(im) .* dV;
        
        % Store data in 'd'
        d.sliceSum(i) = sum(prod1(:));
        d.zVal(i) = Z(i);
        
        % Evaluate time left
        T = toc/60;
        Tremain = round(10*T/i .* (length(sliceList) - i))/10;
        
        % Update status
        disp(['Done ' num2str(i) ' of ' num2str(length(sliceList)) ...
            '   (~' num2str(Tremain) ' min remaining)'])
        
        % Clear variables
        clear im sNum se T Tremain prod1
    end
    
    % Save data
    save([imPath filesep 'I_by_slice'],'d')
    
    % Clear variables
    clear X Y Z Jz Jx Jy dV sliceList
    
end

clear a


%% Calculate I

% Load I data for each slice
load([imPath filesep 'I_by_slice'])
    


function imDisplay(im0)
% Displays image without warnings and with enhanced contrast
warning off
im = adapthisteq(im0);
imshow(im)
warning on
    
function [r,c] = makeMask_sqr(im,idx)
%Used for finding coordinate points on a static image 'bw'.

bw = ones(size(im));

% Create figure window
figure;
set(gcf,'DoubleBuffer','on')
imshow(im)
title('Choose rectangle');
hold on

% Give instructions
disp(' '); disp(' ');
disp('  Choose 2 points to draw a rectangular mask')
disp('  Note: expand figure window for better precision')
disp(' ');
disp('  Left mouse button picks points.');disp(' ');
disp('  Right mouse button removes last point.');disp(' ');
disp('  Press return when done.')

% Initialize variables
n      = 1;
but    = 1;
h      = [];
bwPoly = [];
cMax   = size(im,2);
rMax   = size(im,1);

if idx ==1
    x = 1;
    y = 1;
elseif idx == 2
    x = 1;
    y = rMax;
else
    x = cMax;
    y = rMax;
end

% Interactive loop
while 1 == 1
    [xi,yi,but] = ginput(1);
    
    %Return pressed
    if isempty(but)
        break
        
    %Left click
    elseif but==1
        n = min([2 n+1]);
        x(n) = min([ max([xi 1]) cMax]);
        y(n) = min([ max([yi 1]) rMax]);            
        
    end   
    
    %Update plots
    if ~isempty(h)
        delete(h)
    end
    
    % Define roi
    if idx==1
        [c,r] = meshgrid(1:max(ceil(x)),1:max(ceil(y)));
    elseif idx==2
        [c,r] = meshgrid(1:max(ceil(x)),min(floor(y)):rMax);
    else
        [c,r] = meshgrid(min(floor(x)):cMax,min(floor(y)):rMax);
    end
        
    %c = round([min(x) max(x) max(x) min(x) min(x)]);
    %r = round([max(y) max(y) min(y) min(y) max(y)]);
    
    % Coordinates for plotting
    cP = round([min(x)-1 max(x)+1 max(x)+1 min(x)-1 min(x)-1]);
    rP = round([max(y)+1 max(y)+1 min(y)-1 min(y)-1 max(y)+1]);
    h = fill(cP,rP,[1 0 0]);
    alpha(h,0.5);

end

% Turn r & c into vectors
c = c(:);
r = r(:);

close;