function quick_measure


%% Calibration constants for each of the mag cubes 
% in microns/pixel (using the 63x objective)

cal_1_0X   = 0.2733;
cal_1_6X   = 0.1576;
cal_2_5X   = 0.0984;

units = 'microns';

%% Browse for image file

% Get path
[filename, pathname] = uigetfile({'*'},...
                                  'Select an image file');
% Check input
if isempty(filename)
    return
end

% Load image
im = imread([pathname filesep filename]);

clear filename pathname


%% Prompt for mag cube info

ans = questdlg('Which mag cube?',...
               ' ','1.0','1.6','2.5','2.5');
           
if isempty(ans)
    return
    
elseif strcmp(ans,'1.0')
    calconst = cal_1_0X;
    
elseif strcmp(ans,'1.6')
    calconst = cal_1_6X;
    
elseif strcmp(ans,'2.5')
    calconst = cal_2_5X;
    
end

clear cal_1_0X cal_1_5X cal_2_5X


%% Record points

% Parameters
link = 1;
title_text = 'Record two landmarks';
labels = 0;
numLimit = 2;

% Display image
figure;
if isstruct(im)
    imshow(im.cdata,im.colormap);
else
    imshow(im)
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

len = sqrt(diff(x)^2 + diff(y)^2) * calconst;

disp(' ')
disp(['The length between points is ' num2str(len) ' ' units])
disp(' ')

