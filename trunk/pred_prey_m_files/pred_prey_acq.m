function pred_prey_acq(vPath)

% Acquire kinematics of predator and prey fish


%% Parameter values

% Header for image filenames
nameHead = 'Day';

% Extension for image files
nameSuffix = 'TIF';

% Number of digits for frame number in filename
num_digit = 5;

% Visualize acquisition
visSteps = 0;

% Max number of frames for creating the mean image
maxFrames = 1000;

% Whether to invert the image
invert = 1;


%% Get path of data file, load data

if nargin < 1
    vPath = uigetdir(pwd,'Select first frame');
    if vPath==0
        return
    end
end

% Load filenames for frames
a = dir([vPath filesep  '*' nameHead '*.' nameSuffix]);


%% Define roi

% Look for roi data
a2 = dir([vPath filesep 'roi.mat']);

if isempty(a2)
    
    % Read first frame
    im = imread([vPath filesep a(1).name]);
    
    % Select dimensions of circular roi
    txt = 'Select vertical axis of roi';
    figure;
    [p.roi_v.x,p.roi_v.y]   = choosePoints(im,1,txt);
    
    txt = 'Select horizontal axis of roi';
    [p.roi_h.x,p.roi_h.y]   = choosePoints(im,1,txt);
    
    % Overlay roi
    [x_roi,y_roi] = roiCoords(p);
    hold on
    h = plot(x_roi,y_roi,'r');
    hold off
    
    % Save roi data
    save([vPath filesep 'roi.mat'],'p')
    
    % Clear variables
    clear h x y xVals yVals im a2 roi

end


%% Select starting point

% Look for mean image
a2 = dir([vPath filesep 'start_point.mat']);

if isempty(a2)
    
    % Load p for defined for roi
    load([vPath filesep 'roi.mat'])
    
    warning off all
    
    % Prompt for parameters
    prompt={'Frame rate (per sec)', ...
        'Start frame number',...
        'Last frame number'};
    name='Parameters';
    numlines=1;
    defaultanswer={'10','1',a(end).name(end-num_digit-length(nameSuffix):...
                            end-length(nameSuffix)-1)};
    answer      = inputdlg(prompt,name,numlines,defaultanswer);
    if isempty(answer)
        return
    end
    
    p.framerate = str2num(answer{1});
    startFrame = str2num(answer{2});
    endFrame = str2num(answer{3});
    
    % Get indicies for video frames
    idx = 1;
    for i = 1:length(a)
        frNum = str2num(a(i).name(end-num_digit-length(nameSuffix):...
            end-length(nameSuffix)-1));     
        if (frNum >= startFrame) && (frNum <= endFrame) 
            p.frNums(idx) = frNum;
            p.filename{idx} = a(i).name;
            idx = idx + 1;
        end 
    end
    

    % Clear unneeded
    clear answer prompt name numlines defaultanswer a
    
    % Measure body length, initial position & orientation
    txt = 'Select nose, then caudal peduncle';
    img = imread([vPath filesep p.filename{1}]);
    [xT,yT]   = choosePoints(img,1,txt);
    p.bLength = ((xT(2)-xT(1))^2 + (yT(2)-yT(1))^2)^0.5;
    p.xHead = xT(1);
    p.yHead = yT(1);
    p.xTail = xT(2);
    p.yTail = yT(2);
    p.x = mean(xT);
    p.y = mean(yT);
    clear xT yT txt
    close

    warning on all

    save([vPath filesep 'start_point.mat'],'p');
    
else % if seq_param exists, load

    disp(' '); disp('Loading existing starting point data . . .'); 
    load([vPath filesep 'start_point.mat'])

end

clear img


%% Create or load mean image

% Look for mean image
a2 = dir([vPath filesep 'meanImage.tif']);

% Calculate mean image does not exist
if isempty(a2)   
    
    % Define list of frame numbers, depending on max number of frames
    % requested
    if length(p.frNums) > maxFrames
        dframe = floor(length(p.frNums)/maxFrames);
        frIdx = 1:dframe:length(p.frNums);
        clear dframe
    else
        frIdx = 1:length(p.frNums);
    end
    
    % Create waitbar
    h = waitbar(0,...
            ['Mean image: ' num2str(1)],...
             'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    % Create sum image based on first frame
    [imCurr,tmp] = imread([vPath filesep p.filename{1}]);
    imSum = double(imCurr);
    clear imCurr tmp
    
    
    % Loop through frames 
    for i = 1:length(frIdx)
        
        % Add current frame to sum image
        [imCurr,tmp] = imread([vPath filesep p.filename{frIdx(i)}]);
        imSum        = imSum + double(imCurr);
        clear tmp imCurr
        
        % Update status bar
        h = waitbar(i/length(frIdx),h,...
            ['Mean image: ' num2str(i) ' of ' num2str(length(frIdx)) ' frames']);
        
        % Quit m-file, if cancel button pushed
        if getappdata(h,'canceling')
            close force
            return
        end
        
    end
    
    % Calculate mean from sum image
    imMean = uint8(round(imSum./length(frIdx)));
    
    imMean = imMean(:,:,1);
    
    % Write image to movie dir
    imwrite(imMean,[vPath filesep 'meanImage.tif'],'tif',...
            'Compression','none');
    
    close force
    clear frIdx h i imSum
        
    %imMean = rgb2gray(imMean);
      
    
% Load mean image, if present
else
    
    disp(' ')
    disp('Loading mean image . . .');
    imMean = imread([vPath filesep 'meanImage.tif']);
    
end


%% Select thresholds

% Look for mean image
a2 = dir([vPath filesep 'seq_params.mat']);

if isempty(a2)
    
    % Load p for defined for roi
    load([vPath filesep 'start_point.mat'])
    
    [x_roi,y_roi] = roiCoords(p);
    
    % Grab frames for threshold finding
    img = grabFrame(vPath,p.filename{1},invert,x_roi,y_roi);
    
    % Matlab guesses a threshold value
    p.tVal = graythresh(img);
    
    % Store path info in p
    p.path   = vPath;

    % Run threshFinder to find threshold values
    % note: threshFinder saves p in seq_params.mat
    disp(' ')
    disp('Choose threshold for the predator')
    
    waitfor(threshFinder(img,p))
    

else % if seq_param exists, load

    disp(' '); disp('Loading existing sequence parameters . . .'); 
    load([vPath filesep 'seq_params.mat'])
    disp(' ');

end

clear img


%% Step through frames for position of predator

a3 = dir([vPath filesep 'pred_coords.mat']);

if isempty(a3)
    
    if visSteps
        f = figure;
        set(f,'DoubleBuffer','on')
    end
    
    % Loop through frames
    for i = 1:length(p.frNums)
        
        % Define roi coordinates
        [x_roi,y_roi] = roiCoords(p);
        
        % Grab frame, threshold image & choose roi
        img = grabFrame(vPath,p.filename{i},invert,x_roi,y_roi);
        imROI   = roipoly(img,x_roi,y_roi);
        imBW    = ~im2bw(img,p.tVal);
        imBW    = imBW & imROI;
        clear img
        
        % Get peripheral shapes
        [B,L] = bwboundaries(imBW,'noholes');
        
        % Select blob with greatest periphery
        maxB = 0;
        idx = [];
        for j = 1:length(B)
            if length(B{j}) > maxB
                maxB = length(B{j});
                perim = B{j};
                idx = j;
            end
        end
        
        % Store away data
        pd.frame(i) = p.frNums(i);
        pd.filename{i} = p.filename{i};
        pd.xPerim{i} = perim(:,1);
        pd.yPerim{i} = perim(:,2);
        
        % Visualize frames
        if visSteps
            im = imread([vPath filesep p.filename{i}]);
            figure(f)
            warning off
            imshow(im)
            hold on
            plot(pd.yPerim{end},pd.xPerim{end},'r-')
            title(['Frame ' num2str(p.frNums(i)) ' of ' ...
                num2str(p.frNums(end))])
            hold off
            pause(.2)
            warning on
            clear im
        else
            disp(['Predator acquire: Frame ' num2str(p.frNums(i)) ' of ' ...
                num2str(p.frNums(end))]);
        end
        
        % Clear variables for next loop
        clear img imBW imBW2 props imROI se x_roi y_roi maxB
    end
    
    % Save data
    save([vPath filesep 'pred_coords'],'pd')
    
else
    
    % Load 'pd' structure of predator coordinates
    disp('Loading predator data . . .')
    load([vPath filesep 'pred_coords.mat'])
    
end


%% Acquire position of prey



%TODO: postprocessing on predator to determine head
%TODO: Acquire prey position



function img = grabFrame(dirPath,filename,invert,x_roi,y_roi)
% This version uses a single field from an interlaced video frame
% Note: code could be modified to double temporal resolution

% Load image
img = imread([dirPath filesep filename]);

% Get coordinates for whole frame and individual fields
[X,Y] = meshgrid(1:size(img,2), 1:size(img,1));
[X1,Y1] = meshgrid(1:size(img,2), 1:2:size(img,1));
[X2,Y2] = meshgrid(1:size(img,2), 2:2:size(img,1)-2);

% Extract fields
fld1 = img(1:2:size(img,1),:);
%fld2 = img(2:2:size(img,1),:);

% Interpolate between scan lines of field 1
warning off
fr2 = uint8(interp2(X1,Y1,double(fld1),X2,Y2));
warning on

% Replace field 2 with interpolated values
for i=1:size(fr2,1)
    img(2*i,:) = fr2(i,:);
end

%img = adapthisteq(img,'clipLimit',0.02,'Distribution','rayleigh');

% Load subtraction image
imSub  = imread([dirPath filesep,'meanImage.tif']);   

% Adjust grayscale values and convert to double
im     = (imadjust(img));
imSub  = (imadjust(imSub));

% Subtract background
warning off
im = imsubtract(imSub,im);
warning on

%im(find(im>255))  = 255;

if invert
    im = imcomplement(im);
end

% Use roi to crop image
roiI = roipoly(im,x_roi,y_roi);
img = uint8(255.*ones(size(im,1),size(im,2)));
img(roiI(:)) = im(roiI(:));




function [x_roi,y_roi] = roiCoords(p)
%Provides coordinates for an elliptical region of interest

numPts  = 400;
x_h     = p.roi_h.x(1:2);
y_v     = p.roi_v.y(1:2);
r_h     = abs(x_h(1)-x_h(2))/2;
r_v     = abs(y_v(1)-y_v(2))/2;
x_roi   = [];
y_roi   = [];

theta   = linspace(0,pi/2,round(numPts/4))';
x_roi   = [x_roi; r_h .* cos(theta) + mean(x_h)];
y_roi   = [y_roi; r_v .* sin(theta) + mean(y_v)];

theta   = linspace(pi/2,pi,round(numPts/4))';
x_roi   = [x_roi; r_h .* cos(theta) + mean(x_h)];
y_roi   = [y_roi; r_v .* sin(theta) + mean(y_v)];

theta   = linspace(pi,1.5*pi,round(numPts/4))';
x_roi   = [x_roi; r_h .* cos(theta) + mean(x_h)];
y_roi   = [y_roi; r_v .* sin(theta) + mean(y_v)];

theta   = linspace(1.5*pi,2*pi,round(numPts/4))';
x_roi   = [x_roi; r_h .* cos(theta) + mean(x_h)];
y_roi   = [y_roi; r_v .* sin(theta) + mean(y_v)];


function [x,y] = choosePoints(img,link,txt)
%Used for finding coordinate points on a static image 'img'.
warning off all
imshow(img);
title(txt)
hold on;
set(gcf,'DoubleBuffer','on');
disp(' '); disp(' ');
disp('Left mouse button picks points.');disp(' ');
disp('Right mouse button removes last point.');disp(' ');
disp('Press return to stop.')
n = 0;
but = 1;
while 1 == 1
    [xi,yi,but] = ginput(1);
    if isempty(but)
        break
    elseif but==1
        n = n+1;
        x(n) = xi;
        y(n) = yi;
        if link
            h = plot(x,y,'ro-');
        else
            h = plot(x,y,'ro');
        end
    elseif but==3
        if n-1 < 1
            n = 0;
            x = [];
            y = [];
        else
            n = n-1;
            x = x(1:n);
            y = y(1:n);
        end
        hold off
        imshow(img);
        title(txt)
        hold on
        if link
            h = plot(x,y,'ro-');
        else
            h = plot(x,y,'ro');
        end
    end
end

delete(h)

x = x'; y = y';
warning on all