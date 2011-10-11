function acq_3d(root_path)
% Code used to acquire the 3-d kinematics of larval fish in 3d from 2
% camera views.  This code is designed with impulse chamber experiments in
% mind, but could be adapted to a variety of experimental situations.
% It assumes a right-handed coordinate system, where the z-axis points away
% from gravity and the x-y plane is horizontal.  It also assumes that the
% cameras are directed with perpendicular orientations.
% 
% Each sequence directory should contain the following directories:
%   cal    - should have 2 video frames, each with an image of a
%            ruler/micrometer from the two camera views ('cam1_ruler.tif' 
%            & 'cam2_ruler.tif') and 2 video frames, each including a sharp
%            point that can be viewed from both cams
%            ('cam1_pnt.tif','cam2_pnt.tif')
%   cam1   - video recording from the cam1 view as individual tiffs
%   cam2    - video recording from the cam2 view, as individual tiffs
%



%% General parameters

% Suffix for video frame filenames (e.g. 'tif' or 'TIFF')
im_suf = 'tif';

% Top cam  prefix
pre_camT = 'cam_top_frame_';

% Side cam prefix
pre_camS = 'cam_side_frame_';

% Top cam calibration filename
camT_cal_name = 'cam_top_ruler_001.tif';

% Side cam calibration filename
camS_cal_name = 'cam_side_ruler_001.tif';

% Top cam view of image of same point
camT_pnt_name = 'cam_top_pnt_001.tif';

% Side cam view of image of same point
camS_pnt_name = 'cam_side_pnt_001.tif';


%% Establish directory structure, check all files present

% Prompt for directory
if nargin < 1
   root_path = uigetdir([pwd filesep '*.mat'],'Choose sequence directory');
end

% Define paths for camwera files
camT_path  = [root_path filesep pre_camT];
camS_path  = [root_path filesep pre_camS];

% Paths for calibration images
calT_path = [root_path filesep camT_cal_name];
calS_path = [root_path filesep camS_cal_name];

% Get files for video files
a_camT = dir([camT_path '*.' im_suf]);
a_camS = dir([camS_path '*.' im_suf]);

% Check for "camT" contents
if isempty(a_camT)
    error(['There are no files ending in "' im_suf '" in the directory'])
end

% Check for "camS" contents
if isempty(a_camS)
    error(['There are no files ending in "' im_suf '" in the directory'])
end

% Check for calibration image (top)
if isempty(calT_path)
    error(['No calibration image for ' pre_camT])
end

% Check for calibration image (side)
if isempty(calS_path)
    error(['No calibration image for ' pre_camS])
end


%% Acquire spatial calibration data

if isempty(dir([root_path filesep 'cal_data.mat']))
    
    f = figure;
    
    disp(' ')
    disp('This code assumes a right-handed coordinate system.')
    disp(' The cameras should be arranged as follows:')
    disp(' ')
    disp(' Top cam: y-axis is in the horizontal direction')
    disp('          x-axis is in the vertical direction')
    disp( ' ')
    disp('Side cam: y-axis is in the horizontal direction')
    disp('          z-axis is in the vertical direction')
    disp( ' ')
    
    % Info on planes for top cam ----------------------------------------
    
    % Display video frame for cam1 view
    im_camT = imread([root_path filesep a_camT(1).name]);
    warning off
    imshow(im_camT)
    warning on
    title('Top cam view')
    
    % Prompt for and store direction of y-axis
    title(['Top cam: Click 2 points for the direction of the x-axis'])
    direc = choose_direc(f);
    
    cal.camT.xaxis = direc;
    
    clear direc
    
    % Prompt for and store direction of z-axis
    title(['Top cam: Click 2 points for the direction of the y-axis'])
    direc = choose_direc(f);
    
    cal.camT.yaxis = direc;
    
    clear direc im_camT
    
    
    % Info on planes for side cam ----------------------------------------
    
    % Display video frame for cam2 view
    im_camS = imread([root_path filesep a_camS(1).name]);
    warning off
    imshow(im_camS)
    warning on
    title('Side cam view')
    
    % Prompt for direction of y-axis
    title(['Side cam: Click 2 points for the direction of the y-axis'])
    direc = choose_direc(f);
    
    % Store direction of first axis
    cal.camS.yaxis = direc;
    
    clear direc
    
    % Prompt for direction of z-axis
    title(['Side cam: Click 2 points for the direction of the z-axis'])
    direc = choose_direc(f);
    
    % Store direction of second axis
    cal.camS.zaxis = direc;
    
    clear direc im_calS
    
    close(f)
    
    
    % Calibrate magnificantion in top cam -------------------------------
    im_cal_camT = imread([root_path filesep camT_cal_name]); 
    
    cal.camT.const = calibrate(im_cal_camT,...
            'Top cam: choose 2 points over known disance, press return');
    
    clear im_cal_camT
      
     % Calibrate magnificantion in cam2 -------------------------------  
    im_cal_camS = imread([root_path filesep camS_cal_name]); 
    
    cal.camS.const = calibrate(im_cal_camS,...
              'Side cam: choose 2 points over known disance, press return');
    
    clear im_cal_camS
       
    % Select same point in the two views ------------------------------
    im_pnt_camT = imread([root_path filesep camT_pnt_name]); 
    im_pnt_camS = imread([root_path filesep camS_pnt_name]); 
    
    warning off
    imshow(im_pnt_camT)
    warning on
    
    title('Top cam: Click on point')
    [x,y,b] = ginput(1);
    hold on
    plot(x,y,'r+')
    hold off
    pause(.5)
    
    cal.camT.pnt = [x y];
    
    warning off
    imshow(im_pnt_camS)
    warning on
    
    title('Side cam: Click on point')
    [x,y,b] = ginput(1);
    hold on
    plot(x,y,'r+')
    hold off
    pause(.5)
    
    cal.camS.pnt = [x y];
    
    clear x y b im_pnt_camT im_pnt_camS
    
    
    % Save 'cal'
    save([root_path filesep 'cal_data.mat'],'cal')
    
    clear im_cam1 im_cam2
    close
    
else
    disp(' ')
    disp('Loading calibration data . . .')
    disp(' ')
    
    % Load 'cal' structure
    load([root_path filesep 'cal_data.mat'])
end


%% Acquire points manually

% Create figure window
f = figure;
set(f,'WindowStyle','modal')
set(f,'DoubleBuffer','on')
set(f,'Name','Point selection')

% Read first frame
imT = imread([root_path filesep a_camT(1).name]);
imS = imread([root_path filesep a_camS(1).name]);
    
% Prompt for zoom
warning off
hT(1) = imshow(imT);
title('Top view: select zoom level, press return')
add_labels(cal.camT)
zoom on
pause

% Store top axes in 'd'
d.top_xaxis = get(gca,'XLim');
d.top_yaxis = get(gca,'YLim');

% Prompt for zoom
hS(1) = imshow(imS);
title('Side view: select zoom level, press return')
add_labels(cal.camS)
zoom on
pause

% Store side axes in 'd'
d.side_xaxis = get(gca,'XLim');
d.side_yaxis = get(gca,'YLim');
    
warning on
clear hS hT


% Loop through frames from the two perspectives -------------------------
for i = 1:length(a_camT)
    
    % Read image files for current time 
    imT = imread([root_path filesep a_camT(i).name]);
    imS = imread([root_path filesep a_camS(i).name]);
     
    % Display video frames
    warning off
    hT(1) = imshow(imT);
    title(['Top: click on left eye (frame ' num2str(i) ')'])
    add_labels(cal.camT)
    warning on
    
    % Prompt to input a point 
    [x_tmp,y_tmp,b_tmp] = ginput(1);
    hold on
    plot(x_tmp,y_tmp,'r+')
    hold off
    pause(.2)
    
    % If top view selected
    if get(f,'CurrentObject') == hT(1)
        
        
        
        
        
        %coord(i) = store_coord(coord(i),cal,x_tmp,y_tmp);
        
        hold on 
        hT(2) = plot(x_tmp,y_tmp,'r+');
        
      
    % If side view selected
    else
        
        
    end
    
    
%     if i < length(a_camT)
%         % Store top axes in 'd'
%         d.top_xaxis = set(gca,'XLim');
%         d.top_yaxis = set(gca,'YLim');
% 
%         % Store side axis in 'd'
%         d.side_xaxis = set(gca,'XLim');
%         d.side_yaxis = set(gca,'YLim');
%     end
    
    
    % Reset & delete for next loop
    delete(hT)
    delete(hS)
    
    clear imT imS hT hS
end

        

function frame_to_global(cal,pts,cam_num)
% Transforms video frame coordinates into global coordinate

% Find common axis dimension
if isfield(cal.cam1,'x') && isfield(cal.cam2,'x')
    com_dim = 1;
elseif isfield(cal.cam1,'y') && isfield(cal.cam2,'y')
    com_dim = 2;
elseif isfield(cal.cam1,'z') && isfield(cal.cam2,'z')
    com_dim = 3;
end


function global_to_frame(cal,pts,cam_num)
% Transforms global coordinates into video frame coordinates



function direc = choose_direc(f)
% Calculates a spatial calibration constant (in units/pixel) from user 
% selected points recorded from image of a ruler

figure(f)
set(f,'DoubleBuffer','on');

% Give instructions
disp(' '); disp(' ');
disp('Left mouse button picks points.');disp(' ');
disp('First point selects base.')
disp('Second point selects tip')
disp('')
disp('Press return when done.')

% Prompt for point
[xi,yi,but] = ginput(1);

% Return pressed (quit out)
if isempty(but)
    error('Two points not selected')
else
    x(1) = xi;
    y(1) = yi;
    hold on
    h(1) = plot(x,y,'ro');
end

clear xi yi but

% Prompt for point
[xi,yi,but] = ginput(1);

% Return pressed (quit out)
if isempty(but)
    error('Two points not selected')
end
    
% Determine whether horizontal or vert
if abs(yi-y) > abs(xi-x)
    y(2) = yi;
    x(2) = x(1);
    direc = [0 (y(2)-y(1))/abs(y(2)-y(1))];
elseif abs(yi-y) <= abs(xi-x)
    x(2) = xi;
    y(2) = y(1);
    direc = [(x(2)-x(1))/abs(x(2)-x(1)) 0];
end

% Plot line
h(2) = quiver(x(1),y(1),diff(x),diff(y));
set(h,'Color','r')
pause(.5)
delete(h)
hold off


function S = localSystem(P1,P2,P3)
% Defines a transformation vector for a local coordinate system in an
% inertial frame of reference.  Uses P1 as the xaxis and P2 as the origin, and 
% P3 as the z-axis. Coordinates must be (1x3) vectors. Note: if theses axes 
% are not orthogonal, the z-axis direction is assumed to be more accurate
% than the x-axis and the x-axis direction is adjusted to make the coordinates 
% orthoganal.
 
% Check dimensions of inputs
if size(P1,1)~=1 || size(P1,2)~=3 ||...
   size(P2,1)~=1 || size(P2,2)~=3 ||...
   size(P3,1)~=1 || size(P3,2)~=3
    error('Coordinates must be 1x3 vectors');
end
 
% Define units vectors for x and y axes
xAxis   = (P1-P2)./norm(P1-P2);
zAxis   = (P3-P2)./norm(P3-P2);
 
% Define yaxis from the cross product
yAxis   = cross(zAxis,xAxis);
yAxis   = yAxis./norm(yAxis);
 
% Redefine the xaxis, so all axes are orthoganal
xAxis   = cross(yAxis,zAxis);
 
% Define transformation matrix
S       = [xAxis' yAxis' zAxis'];
 


function [xn,yn,zn] = localToGlobal(x,y,z,origin,S)
% Transforms coordinates from the local coordinate system to the global
% system. Coordinates may be given as nxm matricies of equal dimensions.
 
% Check dimensions of inputs
if ~( (size(origin,1)==1 && size(origin,2)==3) ||...
      (size(origin,1)==3 && size(origin,2)==1) )   
    error('Origin must be a 1x3 or 3x1 vector');
    
elseif size(S,1)~=3 || size(S,2)~=3
    error('S must be a 3x3 matrix');
    
elseif ~min(size(x)==size(y)) || ~min(size(x)==size(z)) || ~min(size(y)==size(z))
    error('x, y, & z must have the same dimensions')
    
end
 
% Loop through column to complete transformation
for i = 1:size(x,2)
    pts     = [x(:,i) y(:,i) z(:,i)];
    pts     = [inv(S)'*pts']';
 
    xn(:,i) = pts(:,1) + origin(1);
    yn(:,i) = pts(:,2) + origin(2);
    zn(:,i) = pts(:,3) + origin(3);
    
    clear pts 
end
 


function [xn,yn,zn] = globalToLocal(x,y,z,origin,S)
% Transforms coordinates from the global coordinate system to the local
% system. Coordinates may be given as nxm matricies of equal dimensions.
 
% Check dimensions of inputs
if ~( (size(origin,1)==1 && size(origin,2)==3) ||...
      (size(origin,1)==3 && size(origin,2)==1) )       
    error('Origin must be a 1x3 or 3x1 vector');
       
elseif size(S,1)~=3 || size(S,2)~=3
    error('S must be a 3x3 matrix');
    
elseif ~min(size(x)==size(y)) || ~min(size(x)==size(z)) || ~min(size(y)==size(z))
    error('x, y, & z must have the same dimensions')
    
end
 
% Loop through column to complete transformation
for i = 1:size(x,2)
    pts         = [x(:,i) y(:,i) z(:,i)];    
    pts(:,1)    = x(:,i)-origin(1);
    pts(:,2)    = y(:,i)-origin(2);
    pts(:,3)    = z(:,i)-origin(3);
    pts         = [S'*pts']';
    
    xn(:,i)     = pts(:,1);
    yn(:,i)     = pts(:,2);
    zn(:,i)     = pts(:,3);
    
    clear pts
end



function add_labels(cam)
% Adds axes to current view, according to calibration data

if isfield(cam,'xaxis')
    if cam.xaxis == [1 0]
        xlabel('x ->')
    elseif cam.xaxis == [-1 0]
        xlabel('<- x')
    elseif cam.xaxis == [0 -1]
        ylabel('x ->')
    elseif cam.xaxis == [0 1]
        ylabel('<- x')
    end
end

if isfield(cam,'yaxis')
    if cam.yaxis == [1 0]
        xlabel('y ->')
    elseif cam.yaxis == [-1 0]
        xlabel('<- y')
    elseif cam.yaxis == [0 -1]
        ylabel('y ->')
    elseif cam.yaxis == [0 1]
        ylabel('<- y')
    end
end

if isfield(cam,'zaxis')
    if cam.zaxis == [1 0]
        xlabel('z ->')
    elseif cam.zaxis == [-1 0]
        xlabel('<- z')
    elseif cam.zaxis == [0 -1]
        ylabel('z ->')
    elseif cam.zaxis == [0 1]
        ylabel('<- z')
    end
end

