function hough_play


%% Parameter values

% Header for image filenames
nameHead = 'seq';

% Extension for image files
nameSuffix = 'tif';

% Visualize acquisition
visSteps = 1;

% Max number of frames for creating the mean image
maxFrames = 1000;

% Whether to invert the image
invert = 1;


%% Get path of data file, load data

if nargin < 1
    vPath = uigetdir(pwd,'Select first frame');
end

% Load filenames for frames
a = dir([vPath filesep  '*' nameHead '*.' nameSuffix]);


%% Define roi

% Look for roi data
a2 = dir([vPath filesep 'roi.mat']);

if isempty(a2)
    
    % Read first frame
    im = imread([vPath filesep a(1).name]);
    
    % Default: image dimensions
    p.roi_v.x = mean(size(im,2)).*[1;1];
    p.roi_v.y = [27; size(im,1)];
    p.roi_h.x = [1;size(im,2)];
    p.roi_h.y = mean(size(im,1)).*[1;1];
    
    % Define roi coordinates
    [x_roi,y_roi] = roiCoords(p);
    
    % Display roi
    figure;
    warning off
    imshow(im)
    warning on
    hold on
    h = plot(x_roi,y_roi,'r');
    hold off
    
    % Prompt for whether the default is good
    answer = questdlg('Will this roi work?','???','Yes','No, define','Cancel');
    
    % Stop, if cancel
    if isempty(answer) || strcmp(answer,'Cancel')
        return
    end
    
    % Proceed, if good
    if strcmp(answer,'Yes')
        % Save roi data
        save([vPath filesep 'roi.mat'],'p')
        
    % Interactively define, if default not good
    elseif strcmp(answer,'No, define')
        
        % Delete prior roi
        delete(h);
        
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
end


%% Create or load mean image

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


%% Select starting point, threshold

% Look for mean image
a2 = dir([vPath filesep 'seq_params.mat']);

if isempty(a2)
    
    % Load p for defined for roi
    load([vPath filesep 'roi.mat'])
    
    yes_okay = 0;

    
    warning off all
    
    % Prompt for parameters
    prompt={'Frame rate (per sec)', ...
        'Start frame number',...
        'Last frame number'};
    name='Parameters';
    numlines=1;
    defaultanswer={'500','1',num2str(length(a))};
    answer      = inputdlg(prompt,name,numlines,defaultanswer);
    if isempty(answer)
        return
    end
    p.framerate = str2num(answer{1});
    p.startFrame = str2num(answer{2});
    p.endFrame = str2num(answer{3});
    
    clear answer prompt name numlines defaultanswer
    
    % Measure body length, initial position & orientation
    txt = 'Select anterior, then posterior ends';
    img = imread([vPath filesep a(p.startFrame).name]);
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
    
    while ~yes_okay
 
        % Grab frames for threshold finding
        img = grabFrame(vPath,a,p.startFrame,invert);
        
        % Matlab guesses a threshold value
        p.tVal = graythresh(img);
        
        % Store path info in p
        p.path   = vPath;
        p.fname  = a(p.startFrame).name;
        
        % Run threshFinder to find threshold values
        % note: threshFinder saves p in seq_params.mat
        disp(' ')
        beep
        disp('Choose threshold for the body (isolated)')
        
        waitfor(threshFinder(img,p))
        load([p.path filesep 'seq_params.mat'])
        
        p.tVal_body = p.tVal;
        
         disp(' ')
        beep
        disp('Choose threshold for the appendages')
        
        waitfor(threshFinder(img,p))
        load([p.path filesep 'seq_params.mat'])
        
        p.tVal_app = p.tVal;
        p = rmfield(p,'tVal');
        
        % Check that it's okay
        [x_roi,y_roi] = roiCoords(p);
        
        figure;
        
        % Image preview
        im      = grabFrame(vPath,a,p.startFrame,invert);
        imROI   = roipoly(im,x_roi,y_roi);
        imBW    = ~im2bw(im,p.tVal_body);
        imBW    = imBW & imROI;
        subplot(1,2,1)
        imshow(imBW)
        title('Sample frame');
        disp('You should be able to see only the larva');  

        ButtonName = questdlg('Is these frames okay?', ...
            'Question', ...
            'Yes - proceed', 'No - redo', 'Cancel', 'Yes - proceed');
        switch ButtonName,
            case 'Yes - proceed',
                yes_okay = 1;
            case 'No - redo',
                yes_okay = 0;
            case 'Cancel',
                return

        end % switch
        
        close
        clear im imROI imBW x_roi y_roi c_on c_off
        warning on all

    end % while loop

else % if seq_param exists, load

    disp(' '); disp('Loading existing sequence parameters . . .'); 
    load([vPath filesep 'seq_params.mat'])

end

clear img


%% Initialize from first frame

roiSize = 2.5;

theta  = -15:15;

fig = figure;


% Load parameter values, p
load([vPath filesep 'seq_params.mat'])

% Create initial image 
im = grabFrame(vPath,a,p.startFrame,invert);

% Body position 
pos = [mean([p.xTail p.xHead]) mean([p.yTail p.yHead])];
 
% Body angle, from coordinates
angl = atan2(p.yHead-pos(2),p.xHead-pos(1));

% Crop and rotate image
[im,roi] = change_image(im,angl,pos,roiSize,p.bLength);

% Calculate displacement
[Dpos,Dangl] = calc_disp(im,theta);

% Store body angle
k.angl = angl + Dangl;

% Store body center
k.xCntr = pos(1)+Dpos*cos(k.angl);
k.yCntr = pos(2)+Dpos*sin(k.angl);

% Store head position
k.xHead = k.xCntr + (p.bLength/2)*cos(k.angl);
k.yHead = k.yCntr + (p.bLength/2)*sin(k.angl);

% Store tail position
k.xTail = k.xCntr - (p.bLength/2)*cos(k.angl);
k.yTail = k.yCntr - (p.bLength/2)*sin(k.angl);



%% Step through frames

frames = p.startFrame:p.endFrame;

for i = 2:length(frames)
   
    % Get video frame
    im = grabFrame(vPath,a,frames(i),invert);

    % Crop and rotate image
%     [im,roi] = change_image(im,k.angl(i-1),[k.xCntr(i-1) k.yCntr(i-1)],...
%                       roiSize,p.bLength);
    [im,roi] = change_image(im,k.angl(i-1),pos,...
                      roiSize,p.bLength);
                  
                  
    % Calculate displacement
    [Dpos,Dangl] = calc_disp(im,theta);
    
    % Store body angle
    k.angl(i,1) = k.angl(i-1,1) + Dangl;
    
    % Store body center
    k.xCntr(i,1) = k.xCntr(i-1)+Dpos*cos(k.angl(i));
    k.yCntr(i,1) = k.yCntr(i-1)+Dpos*sin(k.angl(i));
    
    % Store head position
    k.xHead(i,1) = k.xCntr(i) + (p.bLength/2)*cos(k.angl(i));
    k.yHead(i,1) = k.yCntr(i) + (p.bLength/2)*sin(k.angl(i));
    
    % Store tail position
    k.xTail(i,1) = k.xCntr(i) - (p.bLength/2)*cos(k.angl(i));
    k.yTail(i,1) = k.yCntr(i) - (p.bLength/2)*sin(k.angl(i));

    
    clear Dangl Dpos
    
    
    % Visualize kinematics
    if 0
        warning off
       figure(fig);
       
       subplot(1,2,2)
       imshow(im)
       axis on
       set(gca,'YDir','normal')
       grid on
       
       subplot(1,2,1)
       imshow(grabFrame(vPath,a,frames(i),invert))
       axis on
       set(gca,'YDir','normal')
       hold on
       h = plot(k.xHead(i),k.yHead(i),'ro', ...
           [k.xHead(i) k.xTail(i)],[k.yHead(i) k.yTail(i)],'r-',...
           [roi(1) roi(1)+roi(3) roi(1)+roi(3) roi(1) roi(1)],...
           [roi(2) roi(2) roi(2)+roi(4) roi(2)+roi(4) roi(2)],'k');
       
       hold off
       grid on

       delete(h)
       warning on
    end
    
end

return
    % Visualize transforms
    if 0
 
        figure
        
        subplot(3,2,1)
        imshow(imOld)
        title('Old frame')
        axis on
        
        subplot(3,2,2)
        imshow(imNew)
        title('New frame')
        axis on
        
        subplot(3,2,3)
        imagesc(theta, xpOld, ROld);
        colormap(hot);
        
        xlabel('theta')
        axis square
        
        subplot(3,2,4)
        imagesc(theta, xpNew, RNew);
        colormap(hot);
        xlabel('theta')
        
        axis square
        
        subplot(3,2,5)
        plot(theta,max(ROld))
        xlabel('theta')
        ylabel('R')
        grid on
        
        subplot(3,2,6)
        plot(theta,max(RNew))
        xlabel('theta')
        ylabel('R')
        grid on
        
    end
    
    
    
    
    
    theta = -15:15;
    [R,xp] = radon(imcomplement(imStart),theta);
    
    imagesc(theta, xp, R); colormap(hot);
    xlabel('theta')
    
    % Load image
    %im = rgb2gray(imread([vPath filesep a(i).name]));
   
    
    
    im      = grabFrame(vPath,a,i,invert);
    imBW    = giveBW(im,p);
%     imROI   = roipoly(im,x_roi,y_roi);
%     imBW    = ~im2bw(im,p.tVal);
%     imBW    = imBW & imROI;
    
    stats = regionprops(imBW,im,'Area','Centroid','Orientation',...
                    'WeightedCentroid','MajorAxisLength');
    
    % Get image stats
    %[im,cmap,bw,stats] = visImage(im,imMean,roi,p,0);
    
    % Step through each blob
    if visSteps
        warning off
        imshow(im);
        hold on
        warning on
    end
    
    
    % Find largest blob
    maxArea = 0;
    for j = 1:length(stats)
        if stats(j).Area > maxArea
            iBlob = j;
            maxArea = stats(j).Area;
        end
    end
    
    clear maxArea
    
        
    % Center of area of blob
    xC = stats(iBlob).Centroid(1);
    yC = stats(iBlob).Centroid(2);
    
    % Coordinates of weighted (by pixel intensity) centroid
    xW = stats(iBlob).WeightedCentroid(1);
    yW = stats(iBlob).WeightedCentroid(2);
    
    % Define unit vector for orientation
    unitVect = [xW-xC yW-yC]./norm([xC-xW yC-yW]);
    
    % Approximate the nose position from the major axis length
    xTail = unitVect(1) * 0.9*(stats(iBlob).MajorAxisLength/2) + xW;
    yTail = unitVect(2) * 0.9*(stats(iBlob).MajorAxisLength/2) + yW;
    
    % Plot the results
    if visSteps
        %h = plot(xTail,yTail,'ro');
        %set(h,'MarkerfaceColor','r');
        %hold on
        % = plot([xC xTail],[yC yTail],'r-');
        %set(h,'LineWidth',2)
        %title(['Frame ' num2str(i) ' of ' num2str(length(a))])
    end
    
    % Pause, if visualizing acquisition
    if visSteps
        pause(.05)
        hold off
    else
        disp(['Done ' num2str(i) ' of ' num2str(length(a))])
    end
    



function [im,roi] = change_image(im,angl,pos,roiSize,bLength)
% Crops and rotates images relative to current position and body angle

% Convert angle to degrees
angl = 180*angl/pi;

% Crop image
roi = [pos(1)-roiSize/2*bLength pos(2)-roiSize/2*bLength ...
            roiSize*bLength roiSize*bLength];
        %roi(2)
im  = imcrop(im,roi);

% Rotate image from selected coordinates
% Angle is negative b/c y-axis is flipped
im = imrotate(im,angl,'bilinear','crop');
im(im==0) = 255;

% Mask out edges
mWidth = round(roiSize*bLength/10);
im(1:mWidth,:)         = 255;
im(end-mWidth:end,:)   = 255;
im(:,1:mWidth)         = 255;
im(:,end-mWidth:end)   = 255;


function [Dpos,Dang] = calc_disp(im,theta)
% Performs radon transform to find translational and rotational
% displacement

% Perform radon transform
[R,xp] = radon(imcomplement(im),theta);

% Find rotation from transform
maxR = max(R);
Dang = theta(maxR==max(maxR));

% Find corresponding displacement at that rotation
rVals  = R(:,theta==Dang);
Dpos   = xp(rVals==max(rVals));

% Visualize transform
if 1
    figure
    subplot(2,3,1)
    imagesc(theta, xp, R); colormap(hot);
    xlabel('theta')
    ylabel('xp')
    axis square
    
    subplot(2,3,2:3)
    plot(theta,max(R))
    hold on
    plot([Dang Dang],ylim,'k--')
    xlabel('theta')
    ylabel('max(R)')
    grid on
    
    subplot(2,3,4:5)
    plot(xp,rVals)
    hold on
    plot([Dpos Dpos],ylim,'k--')
    xlabel('xp')
    ylabel('R')
    title(['At theta = ' num2str(Dang)])
    grid on
    xlim([-50 50])
end

% Convert angle to radians
Dang = pi*Dang/180;

if 0
   %subplot(1,2,1)
   imshow(im)
   hold on
   set(gca,'YDir','normal')
   xC = size(im,2)/2;
   yC = size(im,1)/2;
   plot([xC xC+Dpos*cos(Dang)],[yC yC+Dpos*sin(Dang)],'r',...
       [xC+Dpos*cos(Dang)],[yC+Dpos*sin(Dang)],'rx')
   grid on
    axis on
    hold off
    
end





function img = grabFrame(dirPath,a,fNum,invert)

img = imread([dirPath filesep a(fNum).name]);

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

img = im;


% TODO: Add roi




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