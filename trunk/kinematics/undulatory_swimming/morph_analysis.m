function dts = morph_analysis
% Performs a least-squares curve fit for the parabolic shape of the fish's
% head in dorsal and lateral views

%len = quad(@hcurve,0,inf)

im_dir = '/Volumes/Docs/Projects/head_swimming/head_morphometrics/raw_images';
cd(im_dir)
im_files = dir('*.jpg');
load 'im_data.mat'


% Analyzing collected data
j = 1;
for i = 1:7
    dts(i,:) = [i ...
                mean([im_data(j).cal_const*im_data(j).head_length ...
                im_data(j+1).cal_const*im_data(j+1).head_length]) ...
                im_data(j).cal_const*im_data(j).b...
                im_data(j+1).cal_const*im_data(j+1).b];
    
    j = j+2;
end

dts
    
return

%% Run analysis on each sequence

i = length(im_data)+1;
while 1
    
    im = imread(im_files(i).name);
    
    eNum  = 100;

    [x,y]=choosePoints(im,3,0);
    
    ctr = [mean([x(1) x(3)]) mean([y(1) y(3)])];
    
    majoraxis = ((x(2)-ctr(1))^2 + (y(2)-ctr(2))^2)^0.5;
    minoraxis = ((x(1)-ctr(1))^2 + (y(1)-ctr(2))^2)^0.5;
    
    theta = linspace(-pi/2,pi/2,eNum)';
    eX    = majoraxis .* cos(theta)+ctr(1);
    eY    = minoraxis .* sin(theta)+ctr(2);

    imshow(im);
    hold on
    plot(eX,eY,'r')
    title(im_files(i).name)
    hold off
    
    answer = questdlg('Look okay?','Review image','Yes','No','Stop','Yes');
    
    if strcmp(answer,'Yes')
        
        calconst = calibrate(im);
        
        im_data(i).filename    = im_files(i).name;
        im_data(i).coords      = [x y];
        im_data(i).head_length = majoraxis;
        im_data(i).b           = minoraxis;
        im_data(i).ellipse     = [eX eY];
        im_data(i).cal_const   = calconst;
        
        save('im_data','im_data')
        
        i = i+1;
        
    elseif strcmp(answer,'Stop')
        return
    else
        i = i;
    end
    
    close
    
end

%save('im_data','im_data')


function [x,y] = choosePoints(im,numLimit,link,title_text,labels,x_roi,y_roi)
%Used for finding coordinate points on a static image 'im'.
% numLimit is the total number of points desired
% link tells how to diplay the points
% link = 0 - unconnected points
% link = 1 - line drawn between 2 most recent points
% link = 2 - bounding rectangle drawn btwn 2 points
% labels = 1 - adds number labels next to the points
% x_roi & y_roi define a region of interest with 2 points

% Assign defaults for the inputs
if nargin < 5
    labels = 0;
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

warning off

% Check inputs
if link > 2
    link = 0;
    warning('link cannot be greater than 2, set to 0');
end

if (link==2) && ~(numLimit==2)
    numLimit = 2;
    warning('bounding box requires just 2 coordinates, numLimit set to 2');
end

if nargin == 6 
    error('You need to define both x and y coodinates for your roi');
end

if nargin >= 7 && ...
    (~(length(x_roi)==2) || ~(length(y_roi)==2))
    error('Two coordinates are needed to specify a roi');
end
    

% Display image
figure;
if isstruct(im)
    imshow(im.cdata,im.colormap);
else
    imshow(im)
end
if nargin >= 6
    set(gca,'YLim',[min(y_roi) max(y_roi)])
    set(gca,'XLim',[min(x_roi) max(x_roi)])
end
title(title_text)
set(gcf,'DoubleBuffer','on');
 

% Give instructions
disp(' '); disp(' ');
disp('Left mouse button picks points.');disp(' ');
disp('Right mouse button removes last point.');disp(' ');
disp('Press return when done collecting.')
disp('Press esc to exit');
disp(' '); disp(' ');


% Initiate parameter values for loop
n   = 0;
but = 1;
labeloffset = 10;


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
    if nargin >= 6
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
    
    if labels
        for k = 1:length(x)
           ht = text(x(k)+labeloffset,y(k)+labeloffset,num2str(k));
           set(ht,'Color','r')
        end
    end
    
    hold off
    
end
x = x'; y = y';
close;
warning on



function calconst = calibrate(im)
% Calculates a spatial calibration constant (in units/pixel) from user 
% selected points recorded from image of a ruler

% Prompt to load image, if none given
% =========================================================================
if nargin < 1
    [fName, pathName, filterIndex] = uigetfile('*.*', 'Select image file');
    if filterIndex==0
        return
    end
    im = imread([pathName filesep fName]);
end



% Interactively record distance for calibration
% =========================================================================
figure;
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
close

%Calculate distance in pixels
dist_pix = ((x(2)-x(1))^2 + (y(2)-y(1))^2)^0.5;



% Calculate calibration constant
% =========================================================================
answer      = inputdlg('Distance in units:');
dist_units  = str2num(answer{1});

calconst = dist_units ./ dist_pix;

disp(['Calibration constant (units/pix) = ' num2str(calconst)])
