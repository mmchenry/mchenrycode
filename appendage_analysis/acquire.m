function acquire(rootPath)
% Reconstructs the 3d shape of a morphology from 2 side views (assumes
% elliptical cross-sections)
%
% This code requires 5 grayscale images of the same morphology to be saved in the
% same directory.  Each directory should hold the following files:
%   "frontal.jpg" - Image of morph. from top view.
%   "frontal binary.jpg"    - Silhouette of the dorsal view of morph. in white 
%                           (pixval = 255) on an back backdrop (pixval =
%                           0).
%   "side binary.jpg"   - Silhouette of the lateral side of just the body in white (pixval = 255)  
%                           on a black backdrop (pixval = 0).
%   "side.jpg"- Image of the morph. from lateral view.
%
% Code assumes that the side & frontal images have scale bars.



%% Prompt for the filename of images to be analyzed, if not given

if (nargin < 1) 
    rootPath = uigetdir(pwd,'Choose directory for an individual');
end


%% Loading images

imTop = imread([rootPath filesep 'frontal.jpg']);
imLat = imread([rootPath filesep 'side.jpg']);
bwTop = imread([rootPath filesep 'frontal binary.jpg']);
bwLat = imread([rootPath filesep 'side binary.jpg']);


%% Acquiring landmarks

a = dir([rootPath filesep 'data.mat']);

if ~isempty(a)

    disp(' ')
    disp('Loading landmark data (morph) in data.mat . . .')
    load([rootPath filesep 'data.mat'])

else
    
    figure;
    
    % Acquire landmarks from lateral view
    disp(['Trace midline of the morphology using at least 2 points:']);
    disp(['     1. Distal margin at center']);
    disp(['     2 to n-1. Pick off midline points']);
    disp(['     n. Proximal margin at center']);
    
    [morph.lateralLand.x,morph.lateralLand.y]   = choosePoints(bwLat,1);
    
    % Acquire landmarks from top view
    disp(['Choose the same number of points from the top image:']);
    disp(['     1. Distal margin at center']);
    disp(['     2 to n-1. Pick off midline points']);
    disp(['     n. Proximal margin at center']);
    
    [morph.topLand.x,morph.topLand.y]     = choosePoints(bwTop,1);
    
    % Find calibration constants 
    disp(' '); disp('==================================================') 
    disp('Determine calibration constant');
    morph.latCal = calibrate(imLat);
    morph.topCal = calibrate(imTop);
   
    close;
    
    save([rootPath filesep 'data'],'morph')

end
    

%% Finding the peripheral shape of the morphology from the TOP VIEW

% Make image of periphery from bwTop
bwTemp1  = bwperim(im2bw(bwTop,.5));
bwTemp1  = bwTemp1(2:end-1,2:end-1);

% Get coordinates of periphery
[pDorsY,pDorsX]  = find(bwTemp1==1);

% Find rois that separate the two sides of the body
midline.x   = morph.topLand.x';
midline.y   = morph.topLand.y';
imgHeight   = size(bwTop,1);
imgWidth    = size(bwTop,2);

xNodeA      = [1; midline.x; imgWidth; imgWidth; 1];
yNodeA      = [midline.y(1);midline.y; midline.y(end);1; 1];
xNodeB      = [1; midline.x; imgWidth; imgWidth; 1];
yNodeB      = [midline.y(1); midline.y; midline.y(end); imgHeight; imgHeight];

roi1        = roipoly(bwTemp1,xNodeA,yNodeA);
roi2        = roipoly(bwTemp1,xNodeB,yNodeB);

% Get cooridnates of the two peripheral shapes
[morph.periRight.y,morph.periRight.x] = find(  bwTemp1&roi1);
[morph.periLeft.y,morph.periLeft.x]   = find(  bwTemp1&roi2);
  
clear bwTemp1 roi1 roi2 xNodeA yNodeA xNodeB yNodeB imgHeight imgWidth


%% Finding the peripheral shape of the morphology from the LATERAL VIEW

% Make image of periphery from bwLat
bwTemp1     = bwperim(im2bw(bwLat,.5));
bwTemp1     = bwTemp1(2:end-1,2:end-1);
imgHeight   = size(bwLat,1);
imgWidth    = size(bwLat,2);

% Get coordinates of periphery
[pLatY,pLatX] = find(bwTemp1==1);

% Find rois that separate the two sides of the body
midline.x   = morph.lateralLand.x';
midline.y   = morph.lateralLand.y';      

xNodeA      = [1; midline.x; imgWidth; imgWidth; 1];
yNodeA      = [midline.y(1);midline.y; midline.y(end);1; 1];
xNodeB      = [1; midline.x; imgWidth; imgWidth; 1];
yNodeB      = [midline.y(1); midline.y; midline.y(end); imgHeight; imgHeight];

roi1        = roipoly(bwTemp1,xNodeA,yNodeA);
roi2        = roipoly(bwTemp1,xNodeB,yNodeB);

% Get cooridnates of the two peripheral shapes
[morph.periTop.y,morph.periTop.x] = find(  bwTemp1&roi1);
[morph.periBot.y,morph.periBot.x] = find(  bwTemp1&roi2);

clear bwTemp1 roi1 roi2 xNodeA yNodeA xNodeB yNodeB imgHeight imgWidth
  

%% Save and plot data
% This allows visual confirmation that the shape is correctly acquired
  
save([rootPath filesep 'data'],'morph')

warning off

figure;
subplot(2,1,1)
imshow(imLat);
hold on
plot(morph.periTop.x,morph.periTop.y,'g.',...
    morph.periBot.x,morph.periBot.y,'b.')
legend('Top surface','Bottom surface');

subplot(2,1,2)
imshow(imTop);
hold on
plot(morph.periRight.x,morph.periRight.y,'g.',...
    morph.periLeft.x,morph.periLeft.y,'b.')
legend('Right surface','Left surface');

warning on


%% Analyze drag 
% Defines 3rd moment of drag from the acquired data















%% FUNCTIONS

function calconst = calibrate(im)
% Calculates a spatial calibration constant (in units/pixel) from user 
% selected points recorded from image of a ruler

warning off

% Interactively record distance for calibration
imshow(im);

set(gcf,'DoubleBuffer','on');
disp(' '); disp(' ');
disp('Left mouse button picks points.');disp(' ');
disp('Right mouse button removes last point.');disp(' ');
disp('Press return to stop.')

%Loop for interactive input
n = 0; but = 1;
while 1 == 1
    [xi,yi,but] = ginput(1);
    
    %Return pressed (quit out)
    if isempty(but)
        break
        
    %Right click (add point)
    elseif but==1    
        if n + 1 > 2
            n = 2;
        else
            n = n+1;
        end
        
        x(n) = xi;
        y(n) = yi;
        
    %Left click (remove point)
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
    end
    
    % Display points
    imshow(im);
    hold on
    plot(x,y,'ro-'); 
    hold off
end

%Calculate distance in pixels
dist_pix = ((x(2)-x(1))^2 + (y(2)-y(1))^2)^0.5;

% Calculate calibration constant
answer      = inputdlg('Distance in units:');
dist_units  = str2num(answer{1});

calconst = dist_units ./ dist_pix;

warning on


%--------------------------------------------------------------------------
function isThere = isFileThere(fName,pathh)
fName = [fName '.mat'];
isThere = 0;
a = dir(pathh);
for i = 1:length(a)
    tempName    = a(i).name;
    if length(tempName)==length(fName) && ...
            strcmp(tempName,fName)
        isThere = 1;
    end
end


%--------------------------------------------------------------------------   
function [x,y] = choosePoints(img,link)
%Used for finding coordinate points on a static image 'img'.

warning off

imshow(img);
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
            plot(x,y,'ro-')
        else
            plot(x,y,'ro')
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
        hold on
        if link
            plot(x,y,'ro-')
        else
            plot(x,y,'ro')
        end
    end
end
hold off

warning on


%--------------------------------------------------------------------------
function [fileName,filePath,fileExt] = giveFileParts(fileName,pathh)

tmp = find(fileName=='.');

if length(tmp)>1
    error('You cannot have periods in the filename');
elseif length(tmp)<1
    error('Your filename needs an extension')
end

fileExt  = fileName(tmp:end);
fileName = fileName(1:tmp-1);
filePath = [pathh filesep fileName fileExt];


%--------------------------------------------------------------------------
function calConstant = zCalibrateM(im)
f       = figure;
set(f,'Color','k');
set(f,'DoubleBuffer','on');

tVal = 0.5;n=0;i=1;x=[];y=[];
while 1 == 1
    imshow(im);hold on
    h = title(['Frame number = ' num2str(i)]);
    set(h,'Color','w');
    plot(x,y,'r+-')
    [xi,yi,but] = ginput(1);
    if isempty(but)
        break
    elseif but==30
        %increase frame
        if i == numFrames,
            i == numFrames;
        else
            i=i+1;
        end
        x = [];y = [];
    elseif but==46
        %increase frames by 10
        if i+10>numFrames
            i = numFrames;
        else 
            i = i+10;
        end
    elseif but==44
        %decrease frames by 10
        if i-10 <1
            i=1;
        else
            i = i-10;
        end
    elseif but==31
        %decrease frame
        if i==1
            i = 1;
        else
            i = i-1;
        end
        x = [];y = [];
    elseif but==102 | but==70
        %Choose first frame
        i = 1;
    elseif but==109 | but==77
        %Choose middle frame
        i = round(numFrames/2);
    elseif but==108 | but==76
        %Choose middle frame
        i = numFrames;
    elseif but==1
        n = n+1;
        x(n) = xi;
        y(n) = yi;
    elseif but==3
        n = n-1;
        x = x(1:n);
        y = y(1:n);
    end
end
close
if ~length(x)==2
    error(['Incorrect number of points!']);
end
dist        = (diff(x).^2 + diff(y).^2).^.5;
disp(' ');
knownDist   = inputdlg('What is the known distance (note: dont give units)?  ');

calConstant = str2num(knownDist{1}) ./ dist;



