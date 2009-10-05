function calConst = calibrate_image(img)
h           = figure;
            set(h,'Name','Specify distance with two points')
[x,y]       = choosePoints(img,1);
close;
if ~(length(x) == 2)
    error('You need to choose only two points');
end
dist        = inputdlg('What is this distance in SI units? ','Calibration');
calConst    = str2num(dist{1}) ./ norm([x(2)-x(1) y(2)-y(1)]);

function [x,y] = choosePoints(img,link)
%Used for finding coordinate points on a static image 'img'.
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
x = x'; y = y';