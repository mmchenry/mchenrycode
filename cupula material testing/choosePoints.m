function [x,y] = choosePoints(im,numLimit,link,title_text,xlim_c,ylim_c)
%Used for finding coordinate points on a static image 'im'.
% numLimit is the total number of points desired
% link tells how to diplay the points
% link = 0 - unconnected points
% link = 1 - line drawn between 2 most recent points
% link = 2 - bounding rectangle drawn btwn 2 points
% xlim_c,ylim_c define a region of interest with 2 points

% Assign defaults for the inputs
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
    

% Display image
%figure;
if isstruct(im)
    imshow(im.cdata,im.colormap);
else
    imshow(im)
end
if nargin >= 5
    xlim(xlim_c);
    ylim(ylim_c);
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
    if nargin >= 5
        xlim(xlim_c);
        ylim(ylim_c);
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
%close;
