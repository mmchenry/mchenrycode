function measure_bladder
% Interactively prompts user to pick off coordinates and peripheral shape 
% of the swim bladder
 
 
%% Select dimensions from dorsal view (zoomed in)
 
[fName,pName,fIndex] = uigetfile('*.tif','Choose dorsal view, zoomed IN');
if ~fIndex, return; end

cd(pName);
 
im = imread([pName filesep fName]);
[x,y,width] = chooseEllipse(im);
 
if length(x)~=2
    error('You need to select 2 points')
end
 
% Prompt for calConst
format short g;
prompt = {'numerator:','denominator:'};
dlgtitle = 'enter calibration constant Din';

answr = inputdlg(prompt,dlgtitle,1,{'1','1'});
answr = str2double(answr);
answr = answr(1)/answr(2);
%answr = inputdlg({'What is the calibration constant (mm/pix)?'},...
%               'Dorsal, zoomed IN',1,{'.001'});
 

% Stoe data
d.raw.Din_x = x;
d.raw.Din_y = y;
d.raw.Din_width = width;
%d.raw.Din_calConst = str2double(answr{1});
d.raw.Din_calConst = answr;

clear x y width fName fIndex pName answr
 
 
%% Select position from dorsal view (zoomed out)
 
[fName,pName,fIndex] = uigetfile('*.tif','Choose dorsal view, zoomed OUT');
if ~fIndex, return; end
 
im = imread([pName filesep fName]);
[x,y] = choosePoints(im,5,1,1);
 
if length(x)~=5
    error('You need to select 5 points')
end
 
% Prompt for calConst
prompt = {'numerator:','denominator:'};
dlgtitle = 'enter calibration constant Din';

answr = inputdlg(prompt,dlgtitle,1,{'1','1'});
answr = str2double(answr);
answr = answr(1)/answr(2);
%answr = inputdlg({'What is the calibration constant (mm/pix)?'},...
%                'Dorsal, Zoomed out',1,{'.001'});
 
% Store data
d.raw.Dout_x = x;
d.raw.Dout_y = y;
d.raw.Dout_calConst = answr;
 
clear x y fName fIndex pName answ
 
 
%% Select dimensions from lateral view (zoomed in)
 
[fName,pName,fIndex] = uigetfile('*.tif','Choose lateral view, zoomed IN');
if ~fIndex, return; end
 
im = imread([pName filesep fName]);
[x,y,height] = chooseEllipse(im);
 
if length(x)~=2
    error('You need to select 2 points')
end
 
% Prompt for calConst
prompt = {'numerator:','denominator:'};
dlgtitle = 'enter calibration constant Din';

answr = inputdlg(prompt,dlgtitle,1,{'1','1'});
answr = str2double(answr);
answr = answr(1)/answr(2);
%answr = inputdlg({'What is the calibration constant (mm/pix)?'},...
%                'Lateral, Zoomed IN',1,{'.001'});
            
% Store data
d.raw.Lin_x = x;
d.raw.Lin_y = y;
d.raw.Lin_height = height;
d.raw.Lin_calConst = answr;
 
clear x y height fName fIndex answ
 
 
%% Analyze data
 
% Calculate and store position of swim bladder
bodyLength    = ((d.raw.Dout_x(1)-d.raw.Dout_x(5))^2 + ...
                 (d.raw.Dout_y(1)-d.raw.Dout_y(5))^2)^0.5;
bldr_cntr     = [mean([d.raw.Dout_x(2) d.raw.Dout_x(3)]) ...
                 mean([d.raw.Dout_y(2) d.raw.Dout_y(3)])];
nose_ctr      = [d.raw.Dout_x(1) d.raw.Dout_y(1)];
 
d.bldr_position  = ((bldr_cntr(1)-nose_ctr(1))^2 + ...
                    (bldr_cntr(2)-nose_ctr(2))^2)^0.5 / bodyLength;
               
% Calculate volume of swim bladder
a  = 0.5 * ((d.raw.Lin_x(2)-d.raw.Lin_x(1))^2 + ...
            (d.raw.Lin_y(2)-d.raw.Lin_y(1))^2)^0.5 .* d.raw.Lin_calConst; 
c  = d.raw.Din_width.* d.raw.Din_calConst;
b  = d.raw.Lin_height.* d.raw.Lin_calConst;
 
d.bldr_volume = (4/3)*pi*a*b*c;

jj = msgbox(['bladder volume = ' num2str(d.bldr_volume)]);
set(jj,'position',[10 10 150 60]);
 
clear bodyLength bldr_cntr nose_ctr answ
 
 
%% Save data
 
cd(pName);
name = pName(end-10:end);
ii=find(name =='/');
name(ii) = '-';
name = [name '.mat'];
[fName,pName] = uiputfile(name,'Save file');
save([pName filesep fName],'d');

set(jj,'visible','off');
 
 
 
 
 
%% FUNCTION: ChooseEllipse
 
function [x,y,b] = chooseEllipse(im)
%Used for finding an ellipse on a static image 'im'.
 
 
numLimit = 2;
 
warning off
 
% Display image
figure;
 
imshow(im)
set(gcf,'DoubleBuffer','on');
 
 
% Give instructions
disp(' '); disp(' ');
disp('Choose ellipse dimensions in this order:');
disp('1. Posterior point  2. Anterior point   ')
disp('3. Use up/down arrows to select ellipse height');
disp(' ');
disp('Left mouse button picks points.');disp(' ');
disp('Right mouse button removes last point.');disp(' ');
disp('Up/down arrows adjusts height');
disp('Press RETURN when done collecting.')
disp('Press ESC to exit');
disp(' '); disp(' ');
 
 
% Initiate parameter values for loop
n   = 0;
but = 1;
b = [];
bStep = 1;
 
 
% Ellipse coordinates
eNum    = 100;
theta   = linspace(0,2*pi,4*eNum)';
ePts    = [cos(theta) sin(theta)];
 
clear theta 
 
 
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
        
    elseif but==30 % Up arrow
        
        b = b + bStep;
    
    elseif but==31 % Down arrow    
        
        b = b - bStep;
        
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
 
    %title(title_text)
    hold on
    if length(x)==1
        plot(x,y,'+r')
    elseif length(x)==2
        origin = [mean(x) mean(y)];
        S      = localSystem(origin,[x(2) y(2)]);
        a      = 0.5*((x(2)-x(1))^2 + (y(2)-y(1))^2)^0.5;
        bStep  = a/20;
        if isempty(b)
            b = a/2;
        end
        ePtsT  = [a.*ePts(:,1) b.*ePts(:,2)];
        ePtsT  = localToGlobal(ePtsT,origin,S);
        plot(ePtsT(:,1),ePtsT(:,2),'r')
    end
    hold off
    
end
x = x'; y = y';
close;
 
warning on
 
 
%% FUNCTION: ChoosePoints
 
function [x,y] = choosePoints(im,numLimit,link,secondRun)
 
warning off
 
title_text = ' ';
 
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
 
if secondRun
% Give instructions
disp(' '); disp(' ');
disp('Choose points:');
disp('1. Nose   2. Anterior of swim bladder');
disp('3. Posterior of swim bladder   4. Posterior of cellular tail');
disp('5. Posterior of tail test'); 
disp(' ');
disp('Left mouse button picks points.');disp(' ');
disp('Right mouse button removes last point.');disp(' ');
disp('Press return when done collecting.')
disp('Press esc to exit');
disp(' '); disp(' ');
end
 
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
    
    hold off
    
end
x = x'; y = y';
close;
warning on
 
 
 
%% FUNCTIONS: Coordinate systems
 
function S = localSystem(P1,P2)
% Defines a transformation vector for a local coordinate system in an
% inertial frame of reference.  Uses P1 as the origin and P2 to find the
% direction of the x-axis.  Coordinates must be (1x2) vectors.
 
if size(P1,1)~=1 || size(P1,2)~=2 ||...
   size(P2,1)~=1 || size(P2,2)~=2
    error('Coordinates must be (1x2) vectors');
end
 
xAxis       = (P2-P1)./norm(P2-P1);
yAxis       = [xAxis(2);-xAxis(1)];
S           = [xAxis' yAxis];
 
function pts = localToGlobal(pts,origin,S)
 
if size(pts,2)~=2 || size(origin,2)~=2 
    error('Coordinates must be a (nx2) vector');
end
 
pts         = [inv(S)'*pts']';
pts(:,1)    = pts(:,1)+origin(1);
pts(:,2)    = pts(:,2)+origin(2);
 
function pts = globalToLocal(pts,origin,S)
 
if size(pts,2)~=2 || size(origin,2)~=2 
    error('Coordinates must be a (nx2) vector');
end
 
pts(:,1)    = pts(:,1)-origin(1);
pts(:,2)    = pts(:,2)-origin(2);
pts         = [S'*pts']';
