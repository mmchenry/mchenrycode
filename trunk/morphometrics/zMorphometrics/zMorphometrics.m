function zMorphometrics(fileName, batchMode)
% Aquires and analyzes morphology from photographs taken from side and top
% views.  This verison of the code assumes the morphology is the body of a
% larval fish.
%
% batchmode -   (1 or 0) dictates whether to run on all individuals
%
% To conduct the aquistion, 5 grayscale images (alhaving the same filename)
% of the same larva must be saved in the following directories:
%   In "dorsal grayscale"       - Image of the larva from dorsal view.
%   In "dorsal binary"          - Silhouette of the dorsal view of body (no fins)in white 
%                                 (pixval = 255) on an back backdrop (pixval = 0).
%   In "lateral binary_body"    - Silhouette of the lateral side of just the body in white (pixval = 255)  
%                                 on a back backdrop (pixval = 0).
%   in "lateral grayscale"      - Image of the larva from lateral view.
%
% These directories must be saved in zBaseM, defined below
%
% Test message

%% Specify which parts of the code to run

get_raw         = 1;
calc_metrics    = 0;
calc_3d         = 0;
visualize       = 0;

if nargin < 2
    batchMode = 0;
end

%% Directories

zBaseM = '/Users/Bill/Desktop/zMorphometrics_data';
%m_fileDir = [zBaseM filesep 'm_files'];
m_fileDir = '/Users/Bill/Desktop/zMorphometrics';

if (nargin < 1) && ~batchMode
    [fileName,pName,fIndex]= uigetfile({'*.tif';'*.TIF'},'Choose image');
    cd(pName);
    disp(fileName);
    fileName = fileName(1:end-4);
    if ~fIndex
        return
    end
    
    clear fIndex pName
end

%% Parameter values

tolerance   = 1.e4; % Specifies the degree of smoothing of the peripheral shape of the body
numPts_AP   = 200;  % Number of points along the antero-posterior axis
numPts_circ = 200;  % Number of points around the circumference of a cross-section of the body

rho_body    = 1000; % Density of whole body (kg m^-3)
rho_02      = 1;    % Density of oxygen (kg m^-3)
rho_water   = 1000; % Density of water (kg m^-3)
Vbladder    = 1e-10;% Swim bladder volume (in m^3)


%% GET RAW
% Prompts user to calibrate and pick off landmarks from grayscale images, 
% then calculates the peripheral shape of body

if get_raw
    
    % Check batchMode
    if batchMode 
        disp(' '); disp('Running get_raw . . . ');
        files = dir([zBaseM filesep 'to_be_analyzed' filesep ...
                'dorsal grayscale' filesep '*.tif']);
    else
        files(1).name = [fileName '.tif'];
    end
    
    % Loop through files
    for i = 1:length(files)
        
        % Define filename
        fName = files(i).name(1:end-4);
        
        % Load images
        imgDor	  = imread([zBaseM filesep 'to_be_analyzed' filesep ...
            'dorsal grayscale' filesep fName '.tif'],'tif');
        imgLat    = imread([zBaseM filesep 'to_be_analyzed' filesep ...
            'lateral grayscale' filesep fName '.tif'],'tif');
        biDor     = imread([zBaseM filesep 'to_be_analyzed' filesep ...
            'dorsal binary' filesep fName '.tif'],'tif');
        biLat_bod = imread([zBaseM filesep 'to_be_analyzed' filesep ...
            'lateral binary_body' filesep fName '.tif'],'tif');
        
        % Calculate the calibration constant
        tmp = dir([zBaseM filesep 'data_cal' filesep fName '.mat']);
        if isempty(tmp)
            imgCal    = imread([zBaseM filesep 'to_be_analyzed' filesep ...
            'calibration' filesep fName '.tif'],'tif');
            
            
            format short g;
            prompt = {'lat num:','lat den:','dors num:','dors den:','units:'};
            dlgtitle = 'enter calibration constant';

            answer = inputdlg(prompt,dlgtitle,1,{'1','1','1','1','mm'});
            cal.lat.const   = str2num(answer{1}) / str2num(answer{2});
            cal.dors.const  = str2num(answer{3}) / str2num(answer{4});
            cal.units       = answer{5};
            
            % Prompt for calibration constants
            %answer = inputdlg({'Lateral view constant (units/pix)',...
            %                   'Dorsal view constant (units/pix)',...
            %                   'units'},'Calibration constant',...
            %                   1,{'1','1','cm'});
            %cal.lat.const   = str2num(answer{1});
            %cal.dors.const  = str2num(answer{2});
            %cal.units       = answer{3};
            
            % This code interactively determines calbration constant
            %calData = runCalibrations(imgCal,1);
            
            save([zBaseM filesep 'data_cal' filesep ...
                fName],'cal');
            
            clear imgCal calData
        end
        
        clear tmp calData
        
        % If no data, acquire landmarks, find perimeters, separate rois, save data:
        tmp = dir([zBaseM filesep 'data_raw' filesep fileName '.mat']);
        if isempty(tmp)
            
            % Choose landmarks from lateral view
            imSample = imread([m_fileDir filesep 'im_lateral.tiff'],'tif');
            disp(' ');
            messageDisplay('lateral land');
            figure
            [morph.lateralLand.x,morph.lateralLand.y] = ...
                choosePoints(imgLat,1,'Trace midline',imSample);
            close;
            
            if length(morph.lateralLand.x)<2
                error('Choose at least 2 points');
            end
            
            clear imSample
            
            % Choose landmarks from dorsal view
            imSample = imread([m_fileDir filesep 'im_dorsal.tiff'],'tif');
            disp(' ');
            messageDisplay('dorsal land');
            figure;
            [morph.dorsalLand.x,morph.dorsalLand.y] = ...
                choosePoints(imgDor,1,'Trace midline',imSample);
            close;
            
            if length(morph.dorsalLand.x)<2
                error('Choose at least 2 points');
            end
            
            clear imSample
            
            % Choose landmarks of swim bladder from lateral view
            imSample = imread([m_fileDir filesep 'im_swimbladder.tiff'],'tif');
            disp(' ');
            messageDisplay('lateral SB');
            figure;
            [morph.swimBladder.x,morph.swimBladder.y] = ...
                choosePoints(imgLat,1,'Find swim bladder',imSample);
            close;
            
            if length(morph.swimBladder.x)~=2
                error('Choose 2 points for the swim bladder');
            end
            
            clear imSample
            
            % Find coordinates of periphery, using binary images
            morph = givePeriphery(morph,biDor,biLat_bod);
            
            % Save data
            save([zBaseM filesep 'data_raw' filesep fName],'morph');
            
            % Update status
            if batchMode
                disp(' ');
                disp(['    Done ' num2str(i) ' of ' num2str(length(files))])
            end
            
            clear morph
        end
        
        clear tmp fName imgDor imgLat biDor biLat_bod 
    end
end

disp(fileName); %so I know the one I just did
%% CALC METRICS
% Uses the raw data collected to calculate the smoothed shape of the body.
% This will run a batch, if requested

if calc_metrics
    
    % Define files
    if batchMode
        disp(' '); disp('Running calc_metrics . . . ');
        files = dir([zBaseM filesep 'data_raw' filesep '*.mat']);
    else
        files(1).name = fileName;
    end
    
    % Loop through files
    for i = 1:length(files)
        % Define filename
        fName = files(i).name;
        
        % Load data
        load([zBaseM filesep 'data_cal' filesep fName]);
        load([zBaseM filesep 'data_raw' filesep fName]);
        
        % Define calibration constant
        calLat  = cal.lat.const;
        calDors = cal.dors.const;
        units   = cal.units;
        clear cal
        
        % Define periphery of the body (reverse left & right)
        [xR,yR]     = smoothData(morph.periLeft.x,morph.periLeft.y,...
                                    numPts_AP,tolerance);
        [xL,yL]     = smoothData(morph.periRight.x,morph.periRight.y,...
                                    numPts_AP,tolerance);
        [xD,yD]     = smoothData(morph.periDorsal.x,morph.periDorsal.y,...
                                    numPts_AP,tolerance);
        [xV,yV]     = smoothData(morph.periVentral.x,morph.periVentral.y,...
                                    numPts_AP,tolerance);
                                
        % Extract coordinates for the midline (in pix)
        mLat  = [morph.lateralLand.x' morph.lateralLand.y'];
        mDor  = [morph.dorsalLand.x'  morph.dorsalLand.y'];
        
        % Extract position of swim bladder (sb)
        sb_y_pix = mean(morph.swimBladder.x);
        sb_z_pix = mean(morph.swimBladder.y);
        
        % Calculate dimensions of body (SI units)
        w   = abs(yL-yR) .* calDors;
        c   = mean([yD;yV],1) .* calLat;
        h   = abs(yD-yV) .* calLat;
        s   = (xD - xD(1)) .* calDors;
        
        % Calculate position of swim bladder
        sb_y = (sb_y_pix - xL(1)).* calLat;
        sb_z = (sb_z_pix - yL(1)).* calLat;
        sb_x = 0;
        
        % Load grayscale images        
        imDor	  = imread([zBaseM filesep 'to_be_analyzed' filesep ...
            'dorsal grayscale' filesep fName '.tif'],'tif');
        imLat    = imread([zBaseM filesep 'to_be_analyzed' filesep ...
            'lateral grayscale' filesep fName '.tif'],'tif');
        
        clear morph
        
        % Display peripheral shape on grayscale images
        warning off
        figure;
        
        subplot(3,1,1)
        imshow(imDor); hold on
        plot(xL,yL,'r',xR,yR,'b',mDor(:,1),mDor(:,2),'go-')
        legend('left','right')
        
        subplot(3,1,2)
        imshow(imLat); hold on
        plot(xD,yD,'r',xV,yV,'b',sb_y_pix,sb_z_pix,'om',...
             mLat(:,1),mLat(:,2),'go-')
        legend('dorsal','ventral','bladder')
        
        warning on
        
        clear imDor imLat sb_y_pix sb_z_pix
        
        % Plot data
        subplot(3,1,3)
        plot(s,w,'b',s,h,'r',sb_y,sb_z,'ok')
        xlabel(['Body position (' units ')'])
        ylabel(['Length (' units ')'])
        legend('w','h','sb')
                   
        % Calculate trunk step size (check equal steps)
        ds = mean(diff(s));
        if (sum(find(abs(diff(s) - ds) > .001*ds)))
            error('trunk length is not equally spaced');
        end
        
        % Calc segment, total vol of body, vol of tissue
        dV      = pi .* (w./2) .* (h./2) .* ds;
        Vbody   = sum(dV);
        Vtissue = Vbody - Vbladder;  %Vtot was subsituted with Vbody
        
        % Calc mass of the sb, body, tisssue
        Msb     = rho_02 * Vbladder;
        Mbody   = rho_body * Vbody;
        Mtissue = Mbody-Msb;

        % Calculate COV
        COV_z = sum( c.*dV );
        COV_y = sum( s.*dV );
        COV_x = 0;
        COV   = [COV_x COV_y COV_z] ./ Vbody;  %Vbody was Vtot
        
        clear COV_x COV_y COV_z
        
        % Calculate mass/position product for COM (neglecting the bladder)
        COM_z = sum( c.*rho_tissue.*dV );
        COM_y = sum( s.*rho_tissue.*dV );
        COM_x = 0;
        %M     = V .* rho_tissue;
        
        % Subtract tissue within the volume of the swim bladder
        COM_z = COM_z - (rho_tissue*sb_z*vBladder);
        COM_y = COM_y - (rho_tissue*sb_y*vBladder);
        
        % Add mass/position product of swim bladder
        COM_z = COM_z + (sb_z*Msb);
        COM_y = COM_y + (sb_y*Msb);
        
        % Calculate COM from mass/position product
        COM = [COM_x COM_y COM_z]./Mbody;
        
        clear COM_x COM_y COM_z
        
        
        % Store metrics data in m structure
        m.s     = s;        % body position
        m.h     = h;        % height of meat (dist from dorsal to ventral margins)
        m.w     = w;        % width of meat (dist between left and right margins)
        m.c     = c;        % center of meat in verticle dimension
        m.dV    = dV;       % volume at each body segment
        m.Vbody = Vbody;    % total body volume
        m.Vtissue = Vtissue;% Volume of the tissue
        m.Mbody = Mbody;    % Mass of the body
        m.Mtissue = Mtissue;% Mass of the tissue
        m.Msb   = Msb;      % Mass of swim bladder
        m.COV   = COV;      % Center of volume in xyz coordinates
        m.COM   = COM;      % Center of mass in xyz coordinates
        
        % Save
        save([zBaseM filesep 'data_metrics' filesep fName],'m')
        
        % Update status
        if batchMode
            disp(' ');
            disp(['    Done ' num2str(i) ' of ' num2str(length(files))])
        end
        
        clear fName calK s h w c m COV COM dV Vbody Vtissue Mbody Mtissue 
        clear Msb
    end
end

%% CALC 3D
% Uses the raw data collected to calculate the smoothed shape of the body.
% This will run a batch, if requested

if calc_3d
    
    % Define files
    if batchMode
        disp(' '); disp('Running calc_3d . . . ');
        files = dir([zBaseM filesep 'data_metrics' filesep '*.mat']);
    else
        files(1).name = fileName;
    end
    
    % Loop through files
    for i = 1:length(files)
        fName = files(i).name;
        
        % Load data, define dimensions
        load([zBaseM filesep 'data_metrics' filesep fName]);
        
        s = m.s;
        h = m.h;
        w = m.w;
        c = m.c;
        
        clear m
        
        % Define 3d data for visualization
        [s,h,w,c]   = addMouthCap(s,h,w,c);
        [X,Y,Z]     = drawBody(s,h,w,c,numPts_circ);
        
        %Data for visualization:
        threeD.bod.X        = (X-s(1));
        threeD.bod.Y        = Y;
        threeD.bod.Z        = Z;
        
        % Save
        save([zBaseM filesep 'data_3D' filesep fName],'threeD')
        
        % Update status
        if batchMode
            disp(' ');
            disp(['    Done ' num2str(i) ' of ' num2str(length(files))])
        end
        
        clear fName calK s h w c X Y Z threeD
    end
end

%% VISUALIZE
% Visualize the three-dimensional shape of the body:

if visualize
    
    % Check for batch mode
    if batchMode
        disp(' ')
        error('3D visualization only works on individuals, set batchMode = 0');
    end
    
    % Load data
    load([zBaseM filesep 'data_3D' filesep fileName])
    
    % Visualize
    %visualize3D('basic',threeD);
    %visualize3D('for metrics',threeD);
    visualize3D('basic',threeD);
    
end


%% ==========================FUNCTIONS===============================

function morph = givePeriphery(morph,biDor,biLat_bod)
warning off
%Get peripheral points from dorsal view:
  bwTemp1                               = bwperim(im2bw(biDor,.5));
  bwTemp1                               = bwTemp1(2:end-1,2:end-1);
  [pDorsY,pDorsX]                       = find(bwTemp1==1);
  [xNodeA,yNodeA,xNodeB,yNodeB]         = giveNodes(morph,size(biDor,1),size(biDor,2),pDorsX,pDorsY,'dorsal');
  bwTemp2                               = roipoly(im2bw(bwTemp1,.5),xNodeA,yNodeA);
  bwTemp3                               = roipoly(im2bw(bwTemp1,.5),xNodeB,yNodeB);
  [morph.periLeft.y,morph.periLeft.x]   = find(  bwTemp1&bwTemp2);
  [morph.periRight.y,morph.periRight.x] = find(  bwTemp1&bwTemp3);
  
%Get peripheral points of body from lateral view:
  bwTemp1                                   = bwperim(im2bw(biLat_bod,.5));
  bwTemp1                                   = bwTemp1(2:end-1,2:end-1);
  [pLatY,pLatX]                             = find(bwTemp1==1);
  [xNodeA,yNodeA,xNodeB,yNodeB]             = giveNodes(morph,size(biLat_bod,1),size(biLat_bod,2),pLatX,pLatY,'lateral');
  bwTemp2                                   = roipoly(im2bw(bwTemp1,.5),xNodeA,yNodeA);
  bwTemp3                                   = roipoly(im2bw(bwTemp1,.5),xNodeB,yNodeB);
  [morph.periDorsal.y,morph.periDorsal.x]   = find(  bwTemp1&bwTemp2);
  [morph.periVentral.y,morph.periVentral.x] = find(  bwTemp1&bwTemp3);
warning on     
  
function [xNodeA,yNodeA,xNodeB,yNodeB] = giveNodes(morph,imgHeight,imgWidth,pDorsX,pDorsY,action)
% Gives the nodes for a roi that separates the two sides of the body
switch action
case 'dorsal'
    midline.x = morph.dorsalLand.x';
    midline.y = morph.dorsalLand.y';
case 'lateral'
    midline.x = morph.lateralLand.x';
    midline.y = morph.lateralLand.y';      
end
xNodeA      = [1; midline.x; imgWidth; imgWidth; 1];
yNodeA      = [midline.y(1);midline.y; midline.y(end);1; 1];
xNodeB      = [1; midline.x; imgWidth; imgWidth; 1];
yNodeB      = [midline.y(1); midline.y; midline.y(end); imgHeight; imgHeight];
        
function calData = runCalibrations(imgCal,numRuns)
clc;disp('calibrate image');
i = 1;
while i<numRuns+1
    cConstant(i) = zCalibrateM(imgCal);
    if i>1
        if cConstant(i) > 1.1*cConstant(i-1) | cConstant(i) < 0.9*cConstant(i-1)
            clc;disp('WARNING!   2 measurements far apart, starting over');disp(' ');
            i = 0;
        end    
    end
    i = i+1;
end
calData.measurements    = cConstant;
calData.meanCalConst    = mean(cConstant);
calData.precision       = 2.*std(cConstant);      

function messageDisplay(action)
switch action
case 'lateral land'
    disp(['Choose n-landmarks from the image:']);
    disp(['     1. Rostrum']);
    disp(['     2. Center of eye']);
    disp(['     3 to n-2. Pick off midline points']);
    disp(['     n-1. posterior margin of celluar tail']);
    disp(['     n. posterior margin of tail fin']);
    disp(['NOTE!  Be sure to run the midline dorsal to the swim bladder'])
    
case 'lateral SB'
    disp(['Choose n-landmarks from the image:']);
    disp(['     1. Anterior margin of the swim bladder']);
    disp(['     2. Posterior margin of the swim bladder']);
    
case 'dorsal land'
    disp(['Choose n-landmarks from the image:']);
    disp(['     1. Rostrum']);
    disp(['     2. Center between the eyes']);
    disp(['     3 to n-1. Pick off midline points']);
    disp(['     n. ventral margin of celluar tail (ignore fin)']);
end
   
function [x,y] = choosePoints(img,link,title_txt,imSample)
%Used for finding coordinate points on a static image 'img'.

if nargin < 3
    title_txt = ' ';
end

if nargin>3
    subplot(2,1,1)
    imshow(imSample)
    title('Example')
    subplot(2,1,2)
    imshow(img);
    title(title_txt)
else
    imshow(img);
    title(title_txt)
end
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
        if nargin>3
            subplot(2,1,1)
            imshow(imSample)
            title('Example')
            subplot(2,1,2)
            imshow(img);
            title(title_txt)
        else
            imshow(img);
            title(title_txt)
        end
        hold on
        if link
            plot(x,y,'ro-')
        else
            plot(x,y,'ro')
        end
    end
end

function [s,h,w,c] = organizeData(morph,numPts,tolerance)
% Use smoothing spline for each periphery


function [dV,V] = calcVolume(s,h,w,c)
% Works only given a static desciption of body volume
% M - mass of the body
% x,y,z - position of the center of mass in the 'morpho' system

if size(s,1) > size(s,2) | size(h,1) > size(h,2) | size(w,1) > size(w,2)
    error('All inputs must be column vectors');
end


% Calc vertical centroid of area
COMy = 1; 



function [s,h,w,c] = addMouthCap(s,h,w,c)
Ds      = (s(2)-s(1))./10;
s       = [s(1)-2*Ds s(1)-Ds  s];
h       = [h(1)./100 h(1)./2 h];
w       = [w(1)./100 w(1)/2 w];
c       = [c(1) c(1) c];

function [x,y,z]= drawBody(s,h,w,c,numPts)
% Provides 3d coordinates of the body
x=[];y=[];z=[];

theta = linspace(0,2*pi,numPts)';
for i=1:length(s)-1  
  %DRAW FIRST ELLIPSE:
    height      = h(i);
    width       = w(i);
    xTemp1      = s(i)*ones(size(theta));
    yTemp1      = (width/2) .* cos(theta);
    zTemp1      = (height/2) .* sin(theta) + c(i);
    
%     ecc             = ( 1 - (height/2)^2/(width/2)^2 )^.5;
%     [yTemp1,zTemp1] = ellipse1(0,0,[1 ecc],0,[0 360],[],'degrees',numPts);
%   % [yTemp1,zTemp1] = giveEllipse(0,c(i),width/2,height/2,numPts);
%    % [yTemp1,zTemp1] = ellipse1(0,0,[1 axes2ecc(width/2,height/2)],0,[0 360],[],'degrees',numPts);
%     zTemp1          = (width/2) * zTemp1 + c(i);
%     yTemp1          = (width/2) * yTemp1;
    
  %DRAW SECOND ELLIPSE:
    height      = h(i+1);
    width       = w(i+1);
    centerr     = c(i+1);
    
    xTemp2      = s(i+1)*ones(size(theta));
    yTemp2      = (width/2) .* cos(theta);
    zTemp2      = (height/2) .* sin(theta) + c(i+1);
    
%   %  [yTemp2,zTemp2] = giveEllipse(0,c(i+1),width/2,height/2,numPts);
%     ecc             = ( 1 - (height/2)^2/(width/2)^2 )^.5;
%     [yTemp2,zTemp2] = ellipse1(0,0,[1 ecc],0,[0 360],[],'degrees',numPts);
%     %[yTemp2,zTemp2] = ellipse1(0,0,[1 axes2ecc(width/2,height/2)],0,[0 360],[],'degrees',numPts);
%     zTemp2          = (width/2) * zTemp2 + c(i+1);
%     yTemp2          = (width/2) * yTemp2;
%     xTemp2          = s(i+1)*ones(size(yTemp2));   
%     if 0,plot3(xTemp1,yTemp1,zTemp1,'g',xTemp2,yTemp2,zTemp2,'m');axis equal;hold on;end
  %COMBINE DATA:
    all             = 1:length(xTemp2)-1;
    x               = [x [xTemp1(all)'; xTemp2(all)'; xTemp2(all+1)'; xTemp1(all+1)']];
    y               = [y [yTemp1(all)'; yTemp2(all)'; yTemp2(all+1)'; yTemp1(all+1)']];
    z               = [z [zTemp1(all)'; zTemp2(all)'; zTemp2(all+1)'; zTemp1(all+1)']];
end  

function [x,y] = smoothData(x,y,numPts,tolerance)
uniques     = find(~diff(x)==0);
x           = x(uniques);
y           = y(uniques);
sp          = spaps(x, y, tolerance);
x           = min(x):(max(x)-min(x))/(numPts-1):max(x);
y           = fnval(sp,x);
