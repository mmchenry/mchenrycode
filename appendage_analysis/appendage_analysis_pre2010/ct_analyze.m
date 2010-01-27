function ct_analyze(p_name,f_name)
%
% Analyzes ct scans of stomatopod appendages
%
% p_name  - path to ct scna images
% f_name  - filename of first image
%
% Many sections of the code will save results in a mat file in the
% directory containing the images.  These sections of the code will check
% to see if the mat file is already present in the directory.  If so, it
% moves on without executing the code.  If you want it to re-execute, then
% delete the mat file.
%

%% Parameters

num_digits = 3;    % Digits at end of image filenames
scl_factor = 5;    % Factor to downsample image frames (must be integer, >1)
%numslices  = 10;   % Number of slices to make through the propodus
showsteps  = 0;    % Display steps of analysis (1 or 0)
slicesize  = 1;  % Size of slice through propodus (relative to image size)


%% File path

% Prompt for first frame, if not given
if nargin < 2  
    [f_name,p_name,fIndex] = uigetfile({'*.jpg';'*.tif'},...
        'Choose first image in sequence');
    if ~fIndex
        return
    end
end

% Determine number of frames, file information
[tPath,tName,tExt,tVers] = fileparts([p_name filesep f_name]);
prefix = tName(1:end-num_digits);
dirOutput = dir([p_name filesep prefix '*' tExt]);
fileNames = {dirOutput.name}';
numFrames = numel(fileNames);

clear tPath tName tExt tVers dirOutput prefix num_digits


%% Save image data in a structure, mask static text

% Look for data file
tmp = dir([p_name filesep 'ct_data.mat']);

% Execute code, if no data file exists
if isempty(tmp)
      
    % Define sequence of images to import, this translates to a verticle
    % downsampling that is proportional to the degree of xy downsampling 
    im_seq = 1:scl_factor:length(fileNames);

    % Prompt for calibration data    
    prompt={'Image calibration constant (mm/pix):',...
            'Spacing btwn images (mm):'};
    name='Spatial calibration';
    numlines=1;
    defaultanswer={'0.01','0.01'};  
    answer = inputdlg(prompt,name,numlines,defaultanswer);
    
    % Store calibration data
    cal.info       = ['The "full" calibration constants correspond to the '...
                      'full resolution images produced by the scan.  The '...
                      'dwnsmpl values are the calibration constants for '...
                      'the downsampled data.  Downsampling determined by '...
                      'the parameters scl_factor (affects xy data) and '...
                      'img_increment (affect z data).'];
    cal.AR_info    = ['The aspect ratio (AR) is the ratio of the width '...
                      'of the voxel to its depth.'];
    cal.xy_full    = str2num(answer{1});
    cal.AR_full    = (1/cal.xy_full) / (1/str2num(answer{2}));
    cal.xy_dwnsmpl = cal.xy_full*scl_factor;
    
    clear prompt name numlines defaultanswer answer
    
    % The code below assumes AR=1 (if not, need to modify code)
    if cal.AR_full~=1
        error('The resampling in this code assumes AR = 1');
    end
    
    % Read first image, convert to grayscale
    im = imread([p_name filesep fileNames{1}]);
    im = rgb2gray(im);
    im = im(:,:,1);
    
    % Downsample image
    im = imresize(im,1/scl_factor);
        
    % Preallocate the volumetric matrix
    c  = zeros([size(im,1) size(im,2) 1 length(im_seq)],class(im));
    
    % Select masks for static text on frames
    disp(' ')
    disp('Select mask (1 of 3) ====================================')
    [rMsk1,cMsk1] = makeMask_sqr(im);
    
    disp(' ')
    disp('Select mask (2 of 3) ====================================')
    [rMsk2,cMsk2] = makeMask_sqr(im);
    
    disp(' ')
    disp('Select mask (3 of 3) ====================================')
    [rMsk3,cMsk3] = makeMask_sqr(im);
    
    % Check that masks were provided
    if isempty(rMsk1) || isempty(rMsk2) || isempty(rMsk3)
        disp(' '); dips('Quitting: Mask(s) skipped')
        return
    end
            
    j = 1;
    for im_num = im_seq
        
        % Read image, convert to grayscale
        im = imread([p_name filesep fileNames{im_num}]);
        im = rgb2gray(im);
        im = im(:,:,1);

        % Downsample image
        im = imresize(im,1/scl_factor);
        
        % Fill masks
        im(rMsk1,cMsk1) = 0;
        im(rMsk2,cMsk2) = 0;
        im(rMsk3,cMsk3) = 0;
             
        c(:,:,1,j) = im;
        
        j = j+1;
        
        % Update
        disp(['importing image ' num2str(im_num) ' of ' num2str(max(im_seq))])
        
        clear imtmp im_s zT
    end
    
    % Format for isosurface analysis
    disp(' ');
    disp('Calculating smoothed data . . .')
    c = squeeze(c);
    c = smooth3(c);    
    
    % Save data
    save([p_name filesep 'ct_data'],'c');
    save([p_name filesep 'spatial_calibration'],'cal');
    disp(' ');
    disp('Data saved in:')
    disp('ct_data.mat')
    
    clear tmp im_seq c c j im_num cal X Y Z zvals nX nY nZ xT yT
    clear rMsk1 cMsk1 rMsk2 cMsk2 rMsk3 cMsk3 im
end

% Clear variables not needed below
clear scl_factor img_increment f_name numFrames


%% Calibrate pixel values for density

% Look for data file
tmp = dir([p_name filesep 'pixel_calibration.mat']);

% Execute code, if no data file exists
if 0*isempty(tmp)
    
    pix.info = ['This calibration constant is the value you multiply a '...
                'value by to calculate the density (mass/volume) for that '...
                'pixel.'];
            
    % Prompt user for calibration constant and mass of sample
    prompt={'Calibration constant (density/pix val), 0 if unknown:',...
            'Total mass of sample (g):'};
    defaultanswer={'0','0'};  
    answer = inputdlg(prompt,'Pixel value/density calibration',1,defaultanswer);
    
    given_const = str2num(answer{1});
    tot_mass    = str2num(answer{2});
    
    clear prompt defaultanswer answer
    
    % Check answer
    if tot_mass==0
        error('You need to at least enter the measured mass')
    end
    
    % Initialize variable
    cmpr_test = 0;
    
    % Take value of calibration constant, if provided =====================
    if given_const~=0
        pix.cal_const     = given_const;
        pix.mass_measured = tot_mass;
        
        button = questdlg('Compare cal constant calculation (slow)?',...
                          'Validate computations',...
                          'Yes','No','Yes');
        if strcmp(button,'Yes')
            cmpr_test = 1;
        end
        
        clear button given_const
    end
        
    % Otherwise, calculate from the total mass of the item ================
    if given_const==0 || cmpr_test

        % Read first image, convert to grayscale
        im = imread([p_name filesep fileNames{1}]);
        im = rgb2gray(im);
        im = im(:,:,1);

        % Load 'cal', calculate voxel volume
        load([p_name filesep 'spatial_calibration']);
        vox_width = cal.xy_full;
        vox_height= cal.xy_full;
        vox_depth = cal.z_full;
        vox_vol   = vox_width * vox_height * vox_depth;
        pix_area  = vox_width * vox_height;

        % Calculate total volume (in mm^3)
        tot_vol = length(fileNames) * vox_depth * ...
            size(im,1) * vow_width * ...
            size(im,2) * vox_height;
        rho_mean= tot_mass/tot_vol;

        clear vox_width vox_height vox_depth

        % Select masks for static text on frames
        disp(' ')
        disp('Select mask (1 of 3) ====================================')
        [rMsk1,cMsk1] = makeMask_sqr(im);

        disp(' ')
        disp('Select mask (2 of 3) ====================================')
        [rMsk2,cMsk2] = makeMask_sqr(im);

        disp(' ')
        disp('Select mask (3 of 3) ====================================')
        [rMsk3,cMsk3] = makeMask_sqr(im);

        % Check that all masks were provided
        if isempty(rMsk1) || isempty(rMsk2) || isempty(rMsk3)
            disp(' '); dips('Quitting: Mask(s) skipped')
            return
        end

        for i = 1:length(fileNames)

            % Read image, convert to grayscale
            im = imread([p_name filesep fileNames{i}]);
            im = rgb2gray(im);
            im = im(:,:,1);

            % Fill masks
            im(rMsk1,cMsk1) = 0;
            im(rMsk2,cMsk2) = 0;
            im(rMsk3,cMsk3) = 0;

            % Pixel intensity is proportional to (light units)/(pixel area)
            % Therefore, the product of this ratio and the pixel area
            % is equal to the units of light sampled by the pixel, which
            % should be proportional to mass of the voxel.  The sum of 
            % these values is therefore proportional to the total mass 
            % captured by the image, called the light mass.
            lt_mass(i) = sum(im(:).* pix_area);
            
            % Update progress
            disp(['Completed ' num2str(i) ' of ' ...
                  num2str(length(fileNames)) ' images.'])
              
            clear im
        end

        % The mean light density is calculated as the total light mass of
        % all images, divided by the total volume.
        calc_const = rho_mean / sum(lt_mass);

        if cmpr_test
            disp(' ');disp(' ')
            disp('RESULTS OF 2 METHODS FOR CALCULATING CALIBRATION CONSTANT')
            disp('---------------------------------------------------------')
            disp(['Calibration constant provided = ' num2str(pix.cal_const)])
            disp(['Calibration constant calculated = ' num2str(calc_const)])
        end
        
        if given_const==0
            pix.cal_const = calc_const;
        end
        
        clear calc_const rho_mean lt_mass pix_area
        clear rMsk1 cMsk1 rMsk2 cMsk2 rMsk3 cMsk3
    end
    save([p_name filesep 'pixel_calibration'],'pix');
end


%% Interactively select landmarks on the morphology

% Check for presence of data file
tmp = dir([p_name filesep 'landmarks.mat']);

if isempty(tmp)
    disp(' ')
    disp('==========Interactive 3D point selection=============');
    disp(' ')
    disp('Use the Left & Right for the stomatopods body');
    disp(' ')

    % Render 3d image
    load([p_name filesep 'ct_data']);
    h = figure;
    set(gcf,'DoubleBuffer','on');
    h = render3d(c);
    view(3)
    zoom(2)

    description = {'Left on carpus/propodus joint',...
        'Right on carpus/propodus joint',...
        'Distal edge of propodus'};

    for i = 1:length(description)
        disp(' ')
        disp(['Select ' description{i}]);

        [x,y,z] = choosePoints(h,['Pt ' num2str(i)]);

        p(i).coord = [x y z];
        p(i).info = description{i};
        
        clear x y z
    end

    save([p_name filesep 'landmarks'],'p')
    
    clear h i description p
    close
end

clear tmp


%% Transform scan data wrt landmarks

% Check for presence of data file
tmp = dir([p_name filesep 'slice_data.mat']);

if isempty(tmp)

    % Load p and c from data files
    load([p_name filesep 'landmarks.mat'])
    load([p_name filesep 'ct_data.mat']);

    % Define coordinates for c data
    [xc,yc,zc,c] = subvolume(c,[nan,nan,nan,nan,nan,nan]);

    % Define coordinate system from the landmarks
    origin   = mean([p(1).coord; p(2).coord],1);
    z_end_pt = p(3).coord;
    x_end_pt = p(2).coord;
    S        = localSystem(x_end_pt,origin,z_end_pt);

    clear p
    
    % Define number of slices base as the length of the appendage (in
    % pixels)
    numslices = round(norm([z_end_pt-origin]));
    
    % Define a square visual plane for slicing through the ct data
    planesize  = round(slicesize.*max([size(c,1) size(c,2)]));
    z_start    = -norm([z_end_pt-origin])/10;
    z_pos      = linspace(z_start,norm([z_end_pt-origin]),numslices);

    clear slicesize z_start

    % Prep figure if steps to be displayed
    if showsteps
        warning off
        h = figure;
        h = render3d(c,0.5);
        view(121,-16)
        grid on
        
        % Initialize viewing parameters
        tot_views = 10;
        viewIncrement = round(numslices/tot_views);
        viewCount = 0;
    
    end
    
    % Define coordinates for visual slice in local coordinate system
    xPlane  = [1:planesize] - round(planesize/2);
    yPlane  = [1:planesize] - round(planesize/2);
    s.x = xPlane';
    s.y = yPlane';
    s.z = z_pos';

    % Loop through each slice
    for i =1:numslices

        
        zPlane  = zeros(size(xPlane)) + z_pos(i);

        % Format coordinates into matrices
        [xp0,yp0,zp0] = meshgrid(xPlane,yPlane,zPlane);

        % Transform visual plane wrt coordinate system on the propodus
        [xp,yp,zp]    = localToGlobal(xp0,yp0,zp0,origin,S);
        %[xpc,ypc,zpc] = localToGlobal(xpc,ypc,zpc,origin,S);

        % Interpolate for voxel values on visual slice
        ci = interp3(xc,yc,zc,c,xp,yp,zp);

        % Displace visual slices
        if showsteps
            viewCount = viewCount + 1;
            if viewCount==viewIncrement
                %view
                hold on
                h1 = surf(xp,yp,zp,ci);
                set(h1,'LineStyle','none')
                colormap('hot')
                title(['Slice ' num2str(i) ' of ' num2str(numslices)]);
                hold off
                disp(' ')
                disp('Press return to continue')
                pause
                delete(h1)
                clear h1
                viewCount = 0;
            end
        else
            disp(['Slice ' num2str(i) ' of ' num2str(numslices) ' done.'])
        end

        % Store data
        s.info  = 'Coordinates given in the propodus frame of reference';
        s.units = 'pixels';
        s.c(:,:,i)      = ci;
        s.S(:,:,i)      = S;
        s.origin(:,:,i) = origin;

        clear zPlane xpc ypc zpc xp yp zp xpc ypc zpc ci
        clear xp0 yp0 zp0 
    end
    warning on
    close

    save([p_name filesep 'slice_data'],'s')

    clear planesize viewCount tot_views viewIncrement xPlane yPlane
end


%% Select desirable areas from slice data

%tmp = dir([p_name filesep 'slice_data']);

%if isempty(tmp)
    
    numslices = 10;

    %load([p_name filesep 'ct_data']);
    load([p_name filesep 'slice_data'])
    
    % Render surface of scan
    warning off
    hF = figure;
    subplot(1,2,1)
    [X,Y,Z] = meshgrid(s.x,s.y,s.z);
    hF = render3d(s.c,0.5,X,Y,Z);
    view(86,22)
    grid on
    hold on
    
    slices = 1:floor(size(s.c,3)/numslices):size(s.c,3);
    r(1).x = [];r(1).y=[];j = 1;
    
    disp(' ')
    disp(['Note that the number of control points is set by the first '...
         'frame where you specify points'])
    
    while 1
         
        % Collect data for current slice
        curr_slice = slices(j);
        z_pos      = s.z(curr_slice);
        im         = s.c(:,:,curr_slice);
         
        % Format coordinates into matrices
        [xp,yp,zp] = meshgrid(s.x,s.y,z_pos);
        
        
        % Draw slice on 3D image
        subplot(1,2,1)
        h1 = surf(xp,yp,zp,im);
        set(h1,'LineStyle','none')
        colormap('hot')
        %title(['Slice ' num2str(i) ' of ' num2str(numslices)]);
        
        % Select ellipse on slice in 2D
        subplot(1,2,2)
        [r(j).x,r(j).y,cmd] = chooseRoi(im,s.x,s.y,r(j).x,r(j).y,z_pos);
        
        if isempty(cmd) || strcmp(cmd,'next slice')
            j       = j+1;
            
            if j > length(slices)
                break
            else
                r(j).x  = r(j-1).x;
                r(j).y  = r(j-1).y;
            end
            
        elseif strcmp(cmd,'prev slice')
            j = max([1 j-1]);
        
        elseif strcmp(cmd,'reduce slice increment')
            
            % Redefine the spacing in slices
            if j<2
                dslice = round((slices(2)-slices(1))/2);
            else
                dslice = round((slices(j)-slices(j-1))/2);
            end
            slices = [slices(1:j) slices(j)+dslice:dslice:size(s.c,3)];
            
            j       = j+1;
            
            if j > length(slices)
                break
            else
                r(j).x  = r(j-1).x;
                r(j).y  = r(j-1).y;
                
                % Remove any downstream coordinates
                r = r(1:j);
            end
            
       elseif strcmp(cmd,'increase slice increment')
            
            % Redefine the spacing in slices
            if j<2
                dslice = round((slices(2)-slices(1))*2);
            else
                dslice = round((slices(j)-slices(j-1))*2);
            end
            slices = [slices(1:j) slices(j)+dslice:dslice:size(s.c,3)];
            
            j       = j+1;
            
            if j > length(slices)
                break
            else
                r(j).x  = r(j-1).x;
                r(j).y  = r(j-1).y;
                
                % Remove any downstream coordinates
                r = r(1:j);
            end
            
        end
        
        delete(h1)
        
    end
    
    % TODO: create stack that isolates highlighted morphology
    % TODO: Integrate for I calculation
    
%end


%% TODO: Calculate mechanical properties from slice data



%% FUNCTIONS

function [x,y,cmd] = chooseRoi(im,im_x,im_y,x,y,z_pos)
% Used for finding an ellipse on a static image 'im'.
% im_x - x-coodinates of the image
% im_y - y-coordinates of the image
% x & y are the coordinates chosen from a previous trial

if ~isempty(x)
    numlimit = length(x);
else
    numlimit = [];
end

n = length(x);

warning off

% Display image
image(im_x,im_y,im)
axis equal
colormap('gray')
title(['z position = ' num2str(z_pos)]);

set(gcf,'DoubleBuffer','on');
hold on
axis off 

% Give instructions
disp(' '); disp(' ');
disp('Choose region to include with control points');
disp('Start with the center point');
disp(' ');
disp('Control:');
disp('Left mouse - picks points.');disp(' ');
disp('Right mouse - removes last point.');disp(' ');
disp('Up/down arrows - adjust size of ellipse');
disp('Right/Left arrows - move to next/previous slice')
disp('- - Adjust spacing between slices');
disp('e - edits existing points');
disp(' ');
disp('Press RETURN when done collecting.')
disp('Press ESC to exit');
disp(' '); 

% Initiate parameter values for loop

%but = 1;
%b = [];
%bStep = 1;
cmd = [];
a_increment = 0.2;
CCW = 0;
h2 = [];
%h1 = plot(x,y,'r+');
%h2 = plot(x,y,'r+');

editMode = 0;

% Loop for interactive input
while 1 == 1
    
    % Plot ellipse

    h1 = plot(x,y,'or');   
    
    if length(x) > 2
        delete(h2)
        h2 = fill(x,y,'r');
        set(h2,'LineStyle','none')
        alpha(h2,0.5)
    end

    % Get input   
    [xi,yi,but] = ginput(1);
    
    if isempty(but) || but==29
        if length(x)<2
            disp(' ')
            warning(['You only selected ' num2str(length(x)) ' points'])
            
        elseif isempty(numlimit) || length(x)==numlimit  
            break
            
        elseif length(x)<numlimit
            disp(' ')
            disp(['You need to specify ' num2str(numlimit) ' points '...
                     ' before proceeding']);
        end
        
    elseif but==1 % Left click
        
        if editMode
            % Figure which point is closest
            dist = ((xi-x).^2 + (yi-y).^2).^0.5;
            ipnt = find(dist==min(dist),1,'first');

            % Update that point
            x(ipnt,1) = xi;
            y(ipnt,1) = yi;
            
        else
            % If this is run for the first time, don't limit points
            if isempty(numlimit)
                n = n+1;
            
            % Otherwise, limit points
            else
                n = min([numlimit n+1]);
            end
            
            % Record values
            x(n,1) = xi;
            y(n,1) = yi;
            
        end
        
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
        
    elseif but==30 || but==31 % Up arrow/down arrow
        
        if but==30
            chng = 1;
        else
            chng = -1;
        end
        
        if n>2
            cntr = [mean(x) mean(y)];
            displ = chng*a_increment .* [x-cntr(1) y-cntr(2)];
            
            x = x + displ(:,1);
            y = y + displ(:,2);
                
            clear cntr displ
        else
            warning('You need to first define at least 3 points');
        end        
        
    elseif but==101 % The "e" key        
        if ~editMode
            disp(' ')
            disp('EDIT MODE: Choose new location for point')
            disp('EDIT MODE: press "e" again when done')
            editMode = 1;
        else
            disp(' ')
            disp('EDIT MODE COMPLETE')
            editMode = 0;
        end  
        
    elseif but==28 % Left arrow
        
        cmd = 'prev slice';
        break
        
    elseif but==45 % Minus sign
        
        cmd = 'reduce slice increment';
        break
    
    elseif but==45 || but==43 % Plus/equals sign
        
        cmd = 'increase slice increment';
        break
        
    elseif but==27 % If escape
        
        x = [];
        y = [];
        
        break
        
    end

    delete(h1)
    
    
    clear xi yi but
    
    
end

delete(h1)
delete(h2)
hold off
warning on




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coordinate transform functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


function S = localSystem2D(P1,P2)
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


function pts = localToGlobal2D(pts,origin,S)

if size(pts,2)~=2 || size(origin,2)~=2 
    error('Coordinates must be a (nx2) vector');
end

pts         = [inv(S)'*pts']';
pts(:,1)    = pts(:,1)+origin(1);
pts(:,2)    = pts(:,2)+origin(2);


function pts = globalToLocal2D(pts,origin,S)

if size(pts,2)~=2 || size(origin,2)~=2 
    error('Coordinates must be a (nx2) vector');
end

pts(:,1)    = pts(:,1)-origin(1);
pts(:,2)    = pts(:,2)-origin(2);
pts         = [S'*pts']';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions for finding landmarks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x,y,z] = choosePoints(h,title_text)
% Prompts user to select points from 3d surfaces

warning off

% Create empty arrays
x = []; y = []; z =[];

% Set title
title_text_in = title_text;
title([title_text_in ': X = ' num2str(x) ...
       '  Y = ' num2str(y) '  Z = ' num2str(z)]);

while 1

    disp(' ');
    disp('Rotate morphology for best view, press return')
    %figure(h);
    %pause(.01)
    
    hR = rotate3d;
    set(hR,'Enable','on');
    pause
    set(hR,'Enable','off');

    %Select point
    disp(' ');
    disp('Pick point, press return.')
    pause
    p = select3d(gca);
    x = p(1);
    y = p(2);
    z = p(3);

    hold on
    h1 = plot3(x,y,z,'bo');
    set(h1,'MarkerFaceColor','b')
    
    title([title_text_in '  X = ' num2str(x) ...
        '  Y = ' num2str(y) '  Z = ' num2str(z)])
    hold off

    clear p

    disp(' ');
    disp('Rotate morphology to evaluate selection, press return')
    %figure(h);
    
    hR = rotate3d;
    set(hR,'Enable','on');
    pause
    set(hR,'Enable','off');

    button = questdlg('Save point?','Evaluate data',...
        'Yes, proceed','No, try again','Yes, proceed');

    if strcmp(button,'Yes, proceed')
        break
    end
    
    delete(h1)
end

warning on


function h = render3d(c,alphaVal,x,y,z)
% Create surfaces
disp(' ');
disp('Rendering 3D image . . .')

% Find isosurfaces
if nargin > 2
    fvS = isosurface(x,y,z,c);
    fvC = isocaps(x,y,z,c);
    
else
    fvS = isosurface(c);
    fvC = isocaps(c);

end

% Render
hiso = patch(fvS,...
 'FaceColor',.7.*[1 1 1],...
 'EdgeColor','none');

hcap = patch(fvC,...
 'FaceColor','interp',...
 'EdgeColor','none');
colormap('gray')

% Set transparency
if nargin>1
    alpha(hiso,alphaVal);
    alpha(hcap,alphaVal);
end

view(45,30) 
axis equal
%zoom(2)

% Create lights, set renderer, surface properties
l1 = lightangle(45,30); 
l2 = lightangle(-120,0);
set(gcf,'Renderer','OpenGL')
lighting phong
if nargin > 2
    isonormals(x,y,z,c,hiso)
else
    isonormals(c,hiso)
end

set(hcap,'AmbientStrength',.6)
set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50)

xlabel('X');ylabel('Y');zlabel('Z')
title(' ');

h = gcf;


function [r,c] = makeMask_sqr(im)
%Used for finding coordinate points on a static image 'bw'.

bw = ones(size(im));

% Create figure window
figure;
set(gcf,'DoubleBuffer','on')
imshow(im)
title('Choose rectangle');
hold on

% Give instructions
disp(' '); disp(' ');
disp('  Choose 2 points to draw a rectangular mask')
disp('  Note: expand figure window for better precision')
disp(' ');
disp('  Left mouse button picks points.');disp(' ');
disp('  Right mouse button removes last point.');disp(' ');
disp('  Press return when done.')

% Initialize variables
n      = 0;
but    = 1;
h      = [];
bwPoly = [];
cMax   = size(im,2);
rMax   = size(im,1);

% Interactive loop
while 1 == 1
    [xi,yi,but] = ginput(1);
    
    %Return pressed
    if isempty(but)
        break
        
    %Left click
    elseif but==1
        n = min([2 n+1]);
        x(n) = min([ max([xi 1]) cMax]);
        y(n) = min([ max([yi 1]) rMax]);            
        
    %Right click    
    elseif but==3
        n = max([0 n-1]);
        if n==0
            x = [];
            y = [];
        else
            x = x(1:n);
            y = y(1:n);
        end
    end   
    
    %Update plots
    if ~isempty(h)
        delete(h)
    end
    
    % Define roi
    [c,r] = meshgrid(round(min(x):max(x)),round(min(y):max(y)));
    %c = round([min(x) max(x) max(x) min(x) min(x)]);
    %r = round([max(y) max(y) min(y) min(y) max(y)]);
    
    % Coordinates for plotting
    cP = round([min(x)-1 max(x)+1 max(x)+1 min(x)-1 min(x)-1]);
    rP = round([max(y)+1 max(y)+1 min(y)-1 min(y)-1 max(y)+1]);
    h = fill(cP,rP,[1 0 0]);
    alpha(h,0.5);

end

% Turn r & c into vectors
c = c(:);
r = r(:);

close;


function bw = makeMask_poly(im)
%Used for finding coordinate points on a static image 'bw'.

bw = ones(size(im));

lWidth = 2;

subplot(2,1,2)
imshow(bw);
hold on;
subplot(2,1,1)
imshow(im)
title('Make polygon mask');
hold on

set(gcf,'DoubleBuffer','on');

disp(' '); disp(' ');
disp('Choose points to draw a mask on either image')
disp('Note: expand figure window for better precision')
disp(' ');
disp('Left mouse button picks points.');disp(' ');
disp('Right mouse button removes last point.');disp(' ');
disp('Press return when done.')

n = 0;
but = 1;
bwPoly = [];
while 1 == 1
    [xi,yi,but] = ginput(1);
    
    %Return pressed
    if isempty(but)
        break
        
    %Left click
    elseif but==1
        n = n+1;
        x(n) = xi;
        y(n) = yi;
             
        
    %Right click    
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
    
    %Update plots
    subplot(2,1,2)
    bwPoly = roipoly(0,0,im,x,y);
    imshow(~bwPoly & bw)

    subplot(2,1,1)

    imshow(im)
    hold on
    fill(x,y,[1 0 0])
    hold off
    title('Make polygon mask');
end

if ~isempty(bwPoly)
    bw = ~bwPoly & bw;
end
close;

return
% IMAGE PROCESSSING CODE USED IN AN EARLIER VERISON

%         % Surface image: threshold
%         im_s= im2bw(im,.001);
% 
%         % Surface image: fill gaps
%         se = strel('disk',4);
%         im_s = imclose(im_s,se);
%         im_s = imfill(im_s,'holes');
%         
%         % Surface image: find and display boundaries
%         [B,L] = bwboundaries(im_s,'noholes');
%         imtmp = zeros(size(im_s));
%         for k = 1:length(B)
%             bd = B{k};
%             for k2=1:length(bd)
%                 imtmp(bd(k2,1),bd(k2,2)) = 1;
%             end
%            % plot(bd(:,2), bd(:,1), 'r+', 'LineWidth', 1)
%         end