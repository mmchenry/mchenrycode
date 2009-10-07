function zMorphometrics(fileName, batchMode)
% Aquires and analyzes morphology from photographs taken from side and top
% views.  This verison of the code assumes the morphology is the body of a
% larval fish.
%
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


%% Specify which parts of the code to run

get_raw         = 0; % Runs interactive mode for acquiring body periphery
calc_metrics    = 0; % Processes raw data to calc body properties
calc_3d         = 0; % Develop 3d shape of body from data_metrics
visualize3d     = 1; % Visualizes 3d shape
visualizeData   = 0; % Graphs data from data_metrics
simulateFlow    = 0; % Simulates body movement and relative flow

%% Directories

%zBaseM = '/Users/Bill/Desktop/zMorphometrics_data';
zBaseM = '/Volumes/workgroup/relative_velocity/zMorphometrics_data';
%m_fileDir = [zBaseM filesep 'm_files'];
%m_fileDir = '/Users/Bill/Desktop/zMorphometrics';
m_fileDir = '/Volumes/Docs/Projects/Relative velocity/zMorphometrics';

if (nargin < 1) 
    cd(zBaseM)
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

visProfiles = 0;    % Shows traced peripheries on images of larvae

tolerance   = 1.e4; % Specifies the degree of smoothing of the peripheral shape of the body
numPts_AP   = 200;  % Number of points along the antero-posterior axis
numPts_circ = 200;  % Number of points around the circumference of a cross-section of the body

rho_02      = 1;     % Density of oxygen (kg m^-3)
rho_water   = 998;   % Density of water (kg m^-3)  

% Load b, the structure containing body density and swim bladder volume 
load([zBaseM filesep 'body_density_data.mat'])


%% GET RAW
% Prompts user to calibrate and pick off landmarks from grayscale images, 
% then calculates the peripheral shape of body

if get_raw
    
    % Set batchMode
    batchMode = 1;
    
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
    
    disp(fileName); %so I know the one I just did
end


%% CALC METRICS
% Uses the raw data collected to calculate the smoothed shape of the body.
% This will run in batch mode.

if calc_metrics
    
    % Create figure window
    if visProfiles
        hF = figure;
        set(hF,'DoubleBuffer','on')
    end
    
    % Define files
    disp(' '); disp('Running calc_metrics . . . ');
    files = dir([zBaseM filesep 'data_raw' filesep '*.mat']);
    
    % Loop through files
    for i = 1:length(files)
        
        % Define filename, larva number
        fName   = files(i).name(1:end-4);
        larva   = str2num(fName(end-1:end)); 
        
        % Get rho_body (g/mL) and convert to kg m^-3
        rho_body = b.rho_body(larva) * 1000; 
        
        % Get swim bladder volume in mm^3 and convert to m^3 
        Vbladder = b.sb_vol(larva) * 1e-9;   
                                          
        % Load calibration data, in structure 'cal'
        load([zBaseM filesep 'data_cal' filesep fName]);
        
        % Load raw data, stored in structure 'morph'
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
        
        % Calculate dimensions of body (in m)
        w   = abs(yL-yR) .* calDors * 1e-3;
        c   = mean([yD;yV],1) .* calLat * 1e-3;
        h   = abs(yD-yV) .* calLat * 1e-3;
        s   = (xD - xD(1)) .* calDors * 1e-3;
        
        % Extract center position of swim bladder (sb)
        sb_y_pix = mean(morph.swimBladder.x);
        sb_z_pix = mean(morph.swimBladder.y);
        
        % Extract starting and ending positions of swim bladder
        sb_yStart = morph.swimBladder.x(1);
        sb_yEnd   = morph.swimBladder.x(2);
        
        % Calculate position of swim bladder center in m, wrt body
        sb_y = (sb_y_pix - xL(1)).* calLat * 1e-3;
        sb_z = (sb_z_pix - yL(1)).* calLat * 1e-3;
        sb_x = 0;
        
        % Calc starting & ending points on the sb in m, wrt body
        sb_yStart = (sb_yStart - xL(1)).* calLat * 1e-3;
        sb_yEnd   = (sb_yEnd   - xL(1)).* calLat * 1e-3;
               
        % Extract axes for the swim bladder
        sb_B = abs(sb_yEnd - sb_yStart)/2;
        sb_A = sqrt(3*Vbladder/(4*pi*sb_B));
        sb_C = sb_A;
        
        % Load grayscale images        
        imDor	  = imread([zBaseM filesep 'to_be_analyzed' filesep ...
            'dorsal grayscale' filesep fName '.tif'],'tif');
        imLat    = imread([zBaseM filesep 'to_be_analyzed' filesep ...
            'lateral grayscale' filesep fName '.tif'],'tif');
        
        clear morph
        
        % Visualze the results from each individual________________________
        if visProfiles
            
            % Display peripheral shape on grayscale images
            warning off 
            
            subplot(3,1,1)
            imshow(imDor); hold on
            plot(xL,yL,'r',xR,yR,'b',mDor(:,1),mDor(:,2),'go-')
            %legend('left','right')
            title(fName)
            hold off
            
            subplot(3,1,2)
            imshow(imLat); hold on
            plot(xD,yD,'r',xV,yV,'b',sb_y_pix,sb_z_pix,'om',...
                mLat(:,1),mLat(:,2),'go-')
            %legend('dorsal','ventral','bladder')
            hold off
            
            warning on
            
            clear imDor imLat sb_y_pix sb_z_pix
            
            % Plot data
            subplot(3,1,3)
            plot(s,w,'b',s,h,'r',sb_y,sb_z,'ok')
            xlabel(['Body position (' units ')'])
            ylabel(['Length (' units ')'])
            %legend('w','h','sb')
            
            disp('Press return to continue')
            pause
        end
        
        % Calculate trunk step size (check equal steps)
        ds = mean(diff(s));
        if (sum(find(abs(diff(s) - ds) > .001*ds)))
            error('trunk length is not equally spaced');
        end
        
        % Calc seg vol, total vol of body, vol of tissue
        dA      = pi .* (w./2) .* (h./2);
        Vbody   = trapz(s,dA);
        Vtissue = Vbody - Vbladder;  
        
        % Calc mass of the sb, body, tisssue
        Msb        = rho_02 * Vbladder;
        Mbody      = rho_body * Vbody;
        Mtissue    = Mbody-Msb;
        rho_tissue = Mtissue / Vtissue;

        % Calculate COV
        COV_z = trapz(s,c.*dA);
        COV_y = trapz(s,s.*dA);
        COV_x = 0;
        COV   = [COV_x COV_y COV_z] ./ Vbody;  %Vbody was Vtot
        
        clear COV_x COV_y COV_z
        
        
        % COM calculation_____________________________________________
        % Density for just the tissue
        rho_tissue = Mtissue / Vtissue;
        
        % mass-position product for COM (neglecting the bladder)
        COM_z = rho_tissue.*trapz(s,c.*dA);
        COM_y = rho_tissue.*trapz(s,s.*dA);
        COM_x = 0;
        
        % Subtract tissue within the volume of the swim bladder
        COM_z = COM_z - (rho_tissue*sb_z*Vbladder); 
        COM_y = COM_y - (rho_tissue*sb_y*Vbladder);
        
        % Add mass/position product of swim bladder
        COM_z = COM_z + (sb_z*Msb);
        COM_y = COM_y + (sb_y*Msb);
        
        % Calculate COM from mass/position product
        COM = [COM_x COM_y COM_z]./Mbody;
        
        clear COM_x COM_y COM_z

        
        % I calculation_______________________________________________
        % Define positions wrt COM (for I)
        s_COM   = s - COM(2);
        c_COM   = c - COM(3);
        
        % Define values for I analytically integrated wrt x (w/out rho)
        tmp     = 0.25.*(w./2).*(h./2).*pi.*((w./2).^2 + ...
                    4.*(c_COM.^2+s_COM.^2));
        
        % Equation for component of I integrated along x
        % NOTE: tissue in volume of tissue within SB included
        I = rho_tissue.*trapz(s,tmp);
        
        clear tmp s_COM c_COM
        
        % Define vector of transverse sections along y-axis of swim bladder
        s_sb = linspace(sb_yStart,sb_yEnd,500)';
        
        % The axes of elliptical cross-sections along swim bladder
        a0 = real(sqrt(sb_A^2 .* (1 - (s_sb-sb_y).^2./sb_B^2)));
        c0 = real(sqrt(sb_C^2 .* (1 - (s_sb-sb_y).^2./sb_B^2)));
        
        % z-positions of center of dVs in the integration wrt COM
        z0 = sb_z.*ones(size(s_sb)) - COM(3);
        
        % Redefine position of swim bladder sections wrt COM
        s_sb = s_sb - COM(2);
        
        % Equation for component of I integrated along x
        tmp = 0.25.*a0.*c0.*pi.*(a0.^2 + 4.*(z0.^2+s_sb.^2));
        
        % Numerically integrate along sb in y-dimension.  Substract I of
        % tissue & add I of 02 in sb
        I_sb = (rho_02-rho_tissue).*trapz(s_sb,tmp);
        
        clear s_sb a0 c0 z0 tmp
        
        % Substract I of tissue in space occupied by swim bladder & add the
        % I of the air
        I = I + I_sb;
        
        
        % Store metrics data for the individual in m structure_________________
        m.s     = s;        % body position
        m.h     = h;        % height of tissue (dist from dorsal to ventral margins)
        m.w     = w;        % width of tissue (dist between left and right margins)
        m.c     = c;        % center of tissue in verticle dimension
        m.dA    = dA;       % area at each body segment
        m.COV   = COV;      % Center of volume in xyz coordinates
        m.COM   = COM;      % Center of mass in xyz coordinates
        m.I     = I;        % Moment of inertia about COM (z-axis rotation) 
        
        m.info  = {'s: body position' 'h: height of tissue' ...
                   'w: width of tissue' 'c: verticle center of tissue' ...
                   'dA: area at each body segment' ...
                   'COV: Center of volume' ...
                   'COM: Center of mass' ...
                   'I:Moment of inertia about COM (z-axis rotation)'};
                   
        % Store pooled metrics
        mP.b_length(i,1)= s(end);          % Body length
        mP.larvaNum(i,1)= larva;           % Larva number
        mP.age_hr(i,1)  = b.age_hr(larva); % Larva age (hrs)
        mP.rho_body(i,1)= rho_body;        % Body density 
        mP.sb_vol(i,1)  = Vbladder;        % Swim bladder volume
        mP.units        = b.units;
        
        mP.rho_tissue(i,1) = rho_tissue; % Tissue density
        mP.Vbody(i,1)      = Vbody;      % Total body volume
        mP.Vtissue(i,1)    = Vtissue;    % Volume of the tissue
        mP.Mbody(i,1)      = Mbody;      % Mass of the body
        mP.Mtissue(i,1)    = Mtissue;    % Mass of the tissue
        mP.Msb(i,1)        = Msb;        % Mass of swim bladder
        mP.COV(i,:)        = COV;        % Center of volume in xyz coordinates
        mP.COM(i,:)        = COM;        % Center of mass in xyz coordinates
        mP.I(i,1)          = I;          % Moment of inertia about COM (z-axis rotation)
        
        mP.info = {'Vbody: Total vol of body' 'Vtissue: Total vol of tissue' ...
                   'Msb: Mass of swim bladder' 'COV: Center of volume' ...
                   'COM: Center of mass' ...
                   'I:Moment of inertia about COM (z-axis rotation)'};
        
        % Save individual data
        save([zBaseM filesep 'data_metrics' filesep fName],'m')
        
        % Update status
        disp(' ');
        disp(['    Done ' num2str(i) ' of ' num2str(length(files))])
        
        % Clear variables for next pass through the loop
        clear fName calK s h w c m COV COM dV Vbody Vtissue Mbody Mtissue 
        clear Msb larva rho_tissue I 
    end
    
    % Save pooled data
    save([zBaseM filesep 'body_metrics'],'mP')
    
    % Close figure window
    if visProfiles
        close
        clear hF 
    end
end


%% CALC 3D
% Uses the raw data collected to calculate the smoothed shape of the body.
% This will run a batch, if requested

if calc_3d
    
    % Define files
    disp(' '); disp('Running calc_3d . . . ');
    files = dir([zBaseM filesep 'data_metrics' filesep '*.mat']);
    
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
        disp(' ');
        disp(['    Done ' num2str(i) ' of ' num2str(length(files))])

        clear fName calK s h w c X Y Z threeD
    end
end


%% VISUALIZE 3D
% Visualize the three-dimensional shape of the body:

if visualize3d
    
    % Parameters
    alphaVal    = 0.7;
    skinColor   = .5.*[1 1 1];
    backColor   =  1.*[1 1 1];
    COMcolor    = [1 0 0];
    COVcolor    = [0 0 1];
    radiusSize  = .1;        % Radius of COM/COV spheres as a proportion of body width
    sphereRes   = 20;
    
    % Load data threeD data
    load([zBaseM filesep 'data_3D' filesep fileName])
    
    % Load metrics data in structure 'm'
    load([zBaseM filesep 'data_metrics' filesep fileName]);
    
    COV = m.COV;
    COM = m.COM;
    
    clear m
    
    % Create figure
    figure;
    
    % Render body surface
    h1 = patch(real(threeD.bod.X),real(threeD.bod.Y),...
               real(threeD.bod.Z),real(threeD.bod.Z)*0);
    hold on
    set(gcf,'Renderer','OpenGL');
    set(h1,'FaceLighting','gouraud',...
          'LineStyle','none',...
          'BackFaceLighting','reverselit',...
          'FaceColor',skinColor,...
          'AmbientStrength',.5);
    xlabel('x'); ylabel('y');zlabel('z');
    
    % Set transparency of skin
    alpha(h1,alphaVal);
    axis equal;
    
    % Set viewing angle
    viewSide('right');
    
    % Set first light
    l1 = light;
    lightangle(l1,30,60)
    
    % Set second light
    l2 = light;
    lightangle(l2,-30,-60);
    
    % Set figure properties
    set(gcf,'Color',backColor);
    f2= get(gcf,'Children');
    set(f2,'Color',backColor);
    set(f2,'ZDir','reverse');
    
    % Create sphere coordinates
    sphRadius = radiusSize .* range(threeD.bod.X(:))/2;
    [spX,spY,spZ] = sphere(sphereRes);
    
    %Add COM 
    COMx = sphRadius.*spX + COV(1);
    COMy = sphRadius.*spY + COV(2);
    COMz = sphRadius.*spZ + COV(3);
    
    hCOM = surf(COMx,COMy,COMz);
    set(hCOM,'LineStyle','none',...
             'FaceColor',COMcolor);
         
    %Add COV
    COVx = sphRadius.*spX + COM(1);
    COVy = sphRadius.*spY + COM(2);
    COVz = sphRadius.*spZ + COM(3);
    
    hCOV = surf(COVx,COVy,COVz);
    set(hCOV,'LineStyle','none',...
             'FaceColor',COVcolor);
         
    % Add legend
    legend([hCOV;hCOM],{'COV' 'COM'})   
     
    % TODO: Add swim bladder  
    
    % Adjust zoom
    zoom(1)
    %axis off
end


%% Simulate flow
if simulateFlow
    flow_accel = 1.6; % m s^-2
    
    for i = 1:43
        b.rev_vel_20ms(i) = 1000*((1-(0.998/b.rho_body(i)))*1.6*0.02);
        b.rel_vel_20ms_error(i) = 1000*((1-(0.998/(b.rho_body(i)-b.rho_body_error(i))))*1.6*0.02);
        
    end
    
    b.rel_vel_20ms_error = b.rel_vel_20ms - b.rel_vel_20ms_error;

end


%% Visualize data

if visualizeData
    
    % Load pooled data, mP
    load([zBaseM filesep 'body_metrics'])
    age = mP.age_hr./24;
    idx = mP.sb_vol==0;
    
    figure;

    subplot(3,1,1)
    h1 = plot(age(idx),mP.sb_vol(idx),'ko');
    hold on
    set(h1,'MarkerFaceColor','k')
    h2 = plot(age(~idx),mP.sb_vol(~idx),'ko');
    ylabel('Swim bladder volume (m^3)')
    hold off
    
    subplot(3,1,2)
    h1 = plot(age(idx),mP.Vbody(idx),'ko');
    hold on
    set(h1,'MarkerFaceColor','k')
    h2 = plot(age(~idx),mP.Vbody(~idx),'ko');
    ylabel('Body volume (m^3)')
    
    h3 = plot(age(idx),mP.Vtissue(idx),'ro');
    hold on
    %set(h,'MarkerFaceColor','b')
    h4 = plot(age(~idx),mP.Vtissue(~idx),'ro');
    ylabel('Volume (m^3)')    
    hold off
    legend([h2;h3],{'Body' 'Tissue'})
    
    subplot(3,1,3)
    h1 = plot(age(idx),mP.COV(idx,2)./mP.b_length(idx),'ko');
    hold on
    set(h1,'MarkerFaceColor','k')
    h2 = plot(age(~idx),mP.COV(~idx,2)./mP.b_length(~idx),'ko');
    h3 = plot(age(idx),mP.COM(idx,2)./mP.b_length(idx),'r+');
    h4 = plot(age(~idx),mP.COM(~idx,2)./mP.b_length(~idx),'ro');
    
    ylabel('Relative body position')
    legend([h1;h3],{'COV' 'COM'})
    hold off
    
    figure;
    
    subplot(3,1,1)
    h1 = plot(age(idx),mP.Vbody(idx),'ko');
    hold on
    set(h1,'MarkerFaceColor','k')
    h2 = plot(age(~idx),mP.Vbody(~idx),'ko');
    ylabel('Body volume (m^3)')
       
    subplot(3,1,2)
    h1 = plot(age(idx),mP.rho_body(idx)./1000,'ko');
    hold on
    set(h1,'MarkerFaceColor','k')
    h2 = plot(age(~idx),mP.rho_body(~idx)./1000,'ko');
    h3 = plot(age,mP.rho_tissue./1000,'r+');
    legend([h2;h3],{'Body' 'Tissue'})
    ylabel('Density (g/ml)')
    hold off
    
    subplot(3,1,3)
    h1 = plot(age(idx),mP.Mbody(idx),'ko');
    hold on
    set(h1,'MarkerFaceColor','k')
    h2 = plot(age(~idx),mP.Mbody(~idx),'ko');
    ylabel('Body mass (kg)')
    
end



%% FUNCTIONS ======================


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

function viewSide(action)
switch action
case 'dorsal', view([0 90]); case 'ventral', view([0 270]); case 'tail on', view([0 360]);
case 'head on', view([180 360]); case 'right', view([90 0]); case 'left', view([270 0]);
case '3/4 rear', view(3); case '3/4 front',view([40 25]);
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
    yTemp1      = s(i)*ones(size(theta));
    xTemp1      = (width/2) .* cos(theta);
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
    
    yTemp2      = s(i+1)*ones(size(theta));
    xTemp2      = (width/2) .* cos(theta);
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
