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



%TODO: resolve calibration calculation



%% Establish directory structure

% Prompt for directory
if nargin < 1
   root_path = uigetdir([pwd filesep '*.mat'],'Choose sequence directory');
end

% Define cam1 and cam2 paths
cam1_path = [root_path filesep 'cam1']);
cam2_path  = [root_path filesep 'cam2']);

% Check for cam1 cam
if isempty(dir([root_path filesep 'cam1']))
    error(['You need a directory "cam1" that contains your cam1 tiffs'])
end

% Check for cam2 cam
if isempty(dir([root_path filesep 'cam2']))
    error(['You need a directory "cam2" that contains your cam2 tiffs'])
end

% Check for cal
if isempty(dir([root_path filesep 'cal']))
    warning(['There is no "cal" directory that contains calibration tiffs'])
end

% Get files for video files
a_cam1 = dir([cam1_path filesep '*' im_suf]);
a_cam2  = dir([cam2_path filesep '*' im_suf]);

% Check for "cam1" contents
if isempty(a_cam1)
    error('There are no files ending in "' im_suf '" in the cam1 directory'])
end

% Check for "cam2" contents
if isempty(a_cam2)
    error('There are no files ending in "' im_suf '" in the cam2 directory'])
end


%% Spatial calibration

if isempty(dir([root_path filesep 'cal_data.mat']))
    
    figure
    
    
    % Info on planes for cam1 ----------------------------------------
    
    % Display video frame for cam1 view
    im_cam1 = imread([cam1_path filesep a_cam1(1).name]);
    imshow(im_cam1)
    title('cam1 view')
    
    % Gather info on the coordinate system
    button = questdlg('Which plane does cam1 view?','',...
                      'XY','XZ','YZ','XY');
    
    % Store result
    if isempty(button)
        return
    else
        cal.cam1_plane = button;
    end
    
    % Locate origin
    title('Click on the corner of image that is closest to the origin')
    [x,y,b] = ginput(1);
    hold on
    h = plot(x,y,'or');
    hold off
    
    % Determine selected quadrant
    if isempty(b)
        return      
    elseif (x>size(im_cam1,2)/2) && (y<size(im_cam1,1))
        cal.cam1_origin_quad = 1;
    elseif (x<size(im_cam1,2)/2) && (y<size(im_cam1,1))
        cal.cam1_origin_quad = 2;
    elseif (x<size(im_cam1,2)/2) && (y>size(im_cam1,1))
        cal.cam1_origin_quad = 3;
    elseif (x>size(im_cam1,2)/2) && (y>size(im_cam1,1))
        cal.cam1_origin_quad = 4;     
    end
    
    clear x y b h button
    
    
    % Info on planes for cam2 ----------------------------------------
    
    % Display video frame for cam2 view
    im_cam2 = imread([cam2_path filesep a_cam2(1).name]);
    imshow(im_cam2)
    title('cam2 view')
    
    % Gather info on the coordinate system
    button = questdlg('Which plane does cam2 view?','',...
                      'XY','XZ','YZ','XY');
    
    % Store result
    if isempty(button)
        return
    elseif strcmp(cal.cam1_plane,button)
        error('You need cam1 and cam2 to view different planes'])
    else
        cal.cam2_plane = button;
    end
    
    % Locate origin
    title('Click on the corner of image that is closest to the origin')
    [x,y,b] = ginput(1);
    hold on
    h = plot(x,y,'or');
    hold off
    
    % Determine selected quadrant
    if isempty(b)
        return      
    elseif (x>size(im_cam1,2)/2) && (y<size(im_cam1,1))
        cal.cam2_origin_quad = 1;
    elseif (x<size(im_cam1,2)/2) && (y<size(im_cam1,1))
        cal.cam2_origin_quad = 2;
    elseif (x<size(im_cam1,2)/2) && (y>size(im_cam1,1))
        cal.cam2_origin_quad = 3;
    elseif (x>size(im_cam1,2)/2) && (y>size(im_cam1,1))
        cal.cam2_origin_quad = 4;     
    end
    
    clear x y b h button
    
    
    % Confirm that correct orientaion is being calculated ------------
    
    shared_axis = cal.cam2_plane(cal.cam2_plane==cal.cam1_plane);
        
    
    % Calibrate magnificantion in cam1 -------------------------------
    
    
    %cal.cam1 = calibrate
    
    
    
    
    % Calibrate magnificantion in cam2 -------------------------------
    
    clear im_cam1 im_cam2
end




