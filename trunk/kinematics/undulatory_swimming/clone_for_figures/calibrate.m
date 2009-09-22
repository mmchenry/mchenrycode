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