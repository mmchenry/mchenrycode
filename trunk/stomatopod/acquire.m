function acquire(rootPath)
% Reconstructs the 3d shape of a morphology from 2 side views (assumes
% elliptical cross-sections)
%
% This code requires 5 grayscale images of the same morphology to be saved in the
% same directory.  Each directory should hold the following files:
%   "frontal.jpg" - Image of morph. from top view.
%   "frontal binary.jpg"    - Silhouette of the dorsal view of morph. in white 
%                           (pixval = 255) on an back backdrop (pixval =
%                           0).
%   "side binary.jpg"   - Silhouette of the lateral side of just the body in white (pixval = 255)  
%                           on a black backdrop (pixval = 0).
%   "side.jpg"- Image of the morph. from lateral view.
%
% Code assumes that the side & frontal images have scale bars.


%% Prompt for the filename of images to be analyzed, if not given

if (nargin < 1) 
    rootPath = uigetdir(pwd,'Choose directory for an individual');
end


%% Loading images

imTop = imread([rootPath filesep 'frontal.jpg']);
imLat = imread([rootPath filesep 'side.jpg']);
bwTop = imread([rootPath filesep 'frontal binary.jpg']);
bwLat = imread([rootPath filesep 'side binary.jpg']);


%% Acquiring landmarks

a = dir([rootPath filesep 'data_land.mat']);

if ~isempty(a)

    disp(' ')
    disp('Loading landmark data (land) in data.mat . . .')
    load([rootPath filesep 'data_land.mat'])

else
    
    figure;
    
    % Acquire landmarks from lateral view
    clc
    disp('--------------------------------------------------------------');
    disp(['Trace midline of the morphology using at least 2 points:']);
    disp(['     1. Distal margin at center']);
    disp(['     2. Proximal margin at center']);   
    [land.lateralLand.x,land.lateralLand.y]     = choosePoints(bwLat,1);
    
    if length(land.lateralLand.x) ~= 2
        error('Select only 2 points');
    end   
    
    % Acquire axis of rotation
    clc
    disp(' ')
    disp('--------------------------------------------------------------');
    disp(['Select the point of rotation']);   
    [land.lateralRotate.x,land.lateralRotate.y] = choosePoints(imLat,1);
    
    if length(land.lateralRotate.x) ~= 1
        error('Select only 2 points');
    end
    
    % Acquire landmarks from top view
    clc
    disp(['Choose the same number of points from the top image:']);
    disp(['     1. Distal margin at center']);
    disp(['     2. Proximal margin at center']);    
    [land.topLand.x,land.topLand.y]     = choosePoints(bwTop,1);
    
    if length(land.topLand.x) ~= 2
        error('Select only 2 points');
    end 
    
    % Find calibration constants 
    clc
    disp(' '); 
    disp('--------------------------------------------------------------');
    disp('Determine calibration constant');
    land.latCal = calibrate(imLat);
    land.topCal = calibrate(imTop);
   
    close;
    
    save([rootPath filesep 'data_land'],'land')

end


%% Find peripheral coordinates

a = dir([rootPath filesep 'data_edge.mat']);

if ~isempty(a)
    
    disp(' ')
    disp('Loading edge coordinate data in data_edge.mat . . .')
    load([rootPath filesep 'data_edge.mat'])
    
else
    
    numVals  = 100;
    pixRange = 900;
    visData  = 1;
    
    imTop = im2bw(bwTop,.5);
    imLat = im2bw(bwLat,.5);
    
    % Extract data from morph
    latRot  = [land.lateralRotate.x land.lateralRotate.y];
    latLand = [land.lateralLand.x' land.lateralLand.y'];
    topLand = [land.topLand.x' land.topLand.y'];
    
    % Define points for center of elements from lateral view
    xValsLat = linspace(latLand(1,1),latRot(1),numVals)';
    yValsLat = mean(latLand(:,2)).*ones(size(xValsLat,1),size(xValsLat,2));
    
    % Scale length of top view
    sclFctr = abs(diff(topLand(:,1))) / abs(diff(latLand(:,1)));
    
    % Define points for center of elements from top view
    H        = latRot(1) - min(xValsLat);
    xValsTop = sclFctr.*(xValsLat-latLand(1,1)) + topLand(1,1);
    yValsTop = mean(topLand(:,2)).*ones(size(xValsTop,1),size(xValsTop,2));
    
    
    %Set up graph
    if visData
        figure;
        set(gcf,'DoubleBuffer','on')
        %subplot(2,1,1)
        %imshow(imLat)
        %zoom(2)
        %hold on;
        %subplot(2,1,2)
        %imshow(imTop)
        zoom(2)
        hold on
    end
    
    iEnd = find(xValsLat-min(xValsLat) < 0.8.*H,1,'last');
    
    % Loop through points along propodus midline
    for i = 2:iEnd;
        
        
        % LATERAL VIEW ----------------------------------------------------
        
        % Find equation for line that runs perpendicular to vector btwn joint
        % and midline
        slp     = -(xValsLat(i)-latRot(1))/(yValsLat(i)-latRot(2));
        intrcpt = yValsLat(i)-slp.*xValsLat(i);
        
        % Define pixel coordinates along the line
        yPixLat    = linspace(-pixRange/2,pixRange/2,pixRange*2)' + yValsLat(i);
        xPixLat    = (yPixLat-intrcpt) ./ slp;
        
        % Find pixel values along coordinates
        vals       = interp2((imLat),xPixLat,yPixLat);
        cutNum     = round(length(diff(vals))/2);
        
        % Define the first and second half of pixel value vector
        firstHalf  = ([ones(cutNum,1); zeros(length(diff(vals))-cutNum,1)]);
        secondHalf = ([zeros(cutNum,1); ones(length(diff(vals))-cutNum,1)]);
        
        % Find coordinates where there are changes in pixel values
        e1 = [xPixLat(firstHalf & (diff(vals)<0)) ...
            yPixLat(firstHalf & (diff(vals)<0))];
        e2 = [xPixLat(secondHalf & diff(vals)>0) ...
            yPixLat(secondHalf & diff(vals)>0)];
        
        % Store edge coordinates
        edge.frontLat(i-1,:) = e1(1,:);
        edge.backLat(i-1,:)  = e2(1,:);
        edge.midLat(i-1,:)   = [xValsLat(i) yValsLat(i)];
        edge.hVal(i-1,1)     = xValsLat(i) - xValsLat(end);
        edge.rVal(i-1,1)     = ((xValsLat(i)-latRot(1)).^2  +...
                                (yValsLat(i)-latRot(2)).^2)^0.5;
        
        clear vals e1 e2
        
        
        % TOP VIEW --------------------------------------------------------
        
        yPixTop = linspace(-pixRange/2,pixRange/2,pixRange*2)' + yValsTop(i);
        xPixTop = xValsTop(i) .* ones(size(yPixTop,1),size(yPixTop,2));
        
        % Find pixel values along coordinates
        vals       = interp2((imTop),xPixTop,yPixTop);
        
        if min(vals)==max(vals)
            error('No edge detected');
        end
        
        % Find coordinates where there are changes in pixel values
        e1 = [xPixTop(firstHalf & (diff(vals)<0)) ...
            yPixTop(firstHalf & (diff(vals)<0))];
        e2 = [xPixTop(secondHalf & diff(vals)>0) ...
            yPixTop(secondHalf & diff(vals)>0)];
        
        % Store edge coordinates
        edge.frontTop(i-1,:) = e1(1,:);
        edge.backTop(i-1,:)  = e2(1,:);
        edge.midTop(i-1,:)   = [xValsTop(i) yValsTop(i)];
        
        
        
        % Plot data -------------------------------------------------------
        if visData           
            
            %subplot(2,1,1)
            %h1 = plot([xValsLat(i) latRot(1)],[yValsLat(i) latRot(2)],'m-');
            %h2 = plot(xPixLat,yPixLat,'m-');
            %h3 = plot(edge.frontLat(i-1,1),edge.frontLat(i-1,2),'r.');
            %h4 = plot(edge.backLat(i-1,1),edge.backLat(i-1,2),'y.');
            %h5 = plot([edge.frontLat(i-1,1) edge.backLat(i-1,1)],...
            %         [edge.frontLat(i-1,2) edge.backLat(i-1,2)],'r-');
            
            %set(h1,'Color',.2.*[1 1 1])
            %set(h2,'Color',.2.*[1 1 1])
            %clear h1 h2 h3 h4 h5
            
%             subplot(2,1,2)  
             h1 = plot(xPixTop,yPixTop,'b-');
%             h2 = plot(edge.frontTop(i-1,1),edge.frontTop(i-1,2),'r.');
%             h4 = plot(edge.backTop(i-1,1),edge.backTop(i-1,2),'y.');
             h5 = plot([edge.frontTop(i-1,1) edge.backTop(i-1,1)],...
                       [edge.frontTop(i-1,2) edge.backTop(i-1,2)],'r-');
             set(h1,'Color',.2.*[1 1 1])
%             
%             pause(.0001)      
%             clear h1 h2 h3 h4 h5
        end        
        
        disp(['Peripheral coords: done ' num2str(i) ' of ' num2str(iEnd)])
        
    end
    
    if visData
        %subplot(2,1,1)
        plot(latRot(1),latRot(2),'g.')
        hold off
%         subplot(2,1,2)
%         hold off
    end
    
    save([rootPath filesep 'data_edge'],'edge')
    beep; beep
    
end


%% Calculate drag torque index

% Define h, position along dacyl and r, distance to pivot
h = land.latCal .* ...
    [abs(edge.hVal(1))+mean(abs(diff(edge.hVal)));...
     abs(edge.hVal)];
r = land.latCal .* edge.rVal;

% Find values for chord-wise (toward flow) and thickness (perp to flow)
C = land.latCal.* ...
    ((edge.frontLat(:,1) - edge.backLat(:,1)).^2 + ...
     (edge.frontLat(:,2) - edge.backLat(:,2)).^2).^0.5;

T = land.topCal.* ...
    ((edge.frontTop(:,1) - edge.backTop(:,1)).^2 + ...
     (edge.frontTop(:,2) - edge.backTop(:,2)).^2).^0.5;
  
% Calculate drag coefficient
Cd = 0.015.*(1+C./T) + 1.1.*(T./C);

% Sort ascending order of data
idx = length(h)-1:-1:1;
h   = h(length(h):-1:1);
Cd  = Cd(idx);
T   = T(idx);
C   = C(idx);
r   = r(idx);

% Calculate the torque drag index
D = trapz(h.*1000, [Cd.*(T.*1000).*(r.*1000).^3; 0])./((1000.*range(h)).^5);

disp(['D (mm^5) = ' num2str(D)])

% Plot
figure;
subplot(3,1,1)
    plot(h(1:end-1).*1000,C.*1000)
    hold on
    plot(h(1:end-1).*1000,T.*1000,'r')
    legend('C','T')
    ylabel('length (mm)')
subplot(3,1,2)
    plot(h(1:end-1).*1000,Cd)
    ylabel('Cd')
subplot(3,1,3)
    plot(h.*1000,[Cd.*(T.*1000).*(r.*1000).^3; 0],'-')
    xlabel('h (mm)');
    ylabel('Cd T r^3 (mm^5)')
    title(['D = ' num2str(D)]);

% Store and save data
drg.h  = h;
drg.C  = C;
drg.T  = T;
drg.Cd = Cd;
drg.r  = r;
drg.D  = D;
    
save([rootPath filesep 'drag_morph'],'drg')
    

return


%% Save and plot data
% This allows visual confirmation that the shape is correctly acquired
  
save([rootPath filesep 'data'],'morph')

warning off

figure;
subplot(2,1,1)
imshow(imLat);
hold on
plot(land.periTop.x,land.periTop.y,'g.',...
    land.periBot.x,land.periBot.y,'b.')
legend('Top surface','Bottom surface');

subplot(2,1,2)
imshow(imTop);
hold on
plot(land.periRight.x,land.periRight.y,'g.',...
    land.periLeft.x,land.periLeft.y,'b.')
legend('Right surface','Left surface');

warning on


%% Analyze drag 
% Defines 3rd moment of drag from the acquired data















%% FUNCTIONS

function calconst = calibrate(im)
% Calculates a spatial calibration constant (in units/pixel) from user 
% selected points recorded from image of a ruler

warning off

% Interactively record distance for calibration
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

%Calculate distance in pixels
dist_pix = ((x(2)-x(1))^2 + (y(2)-y(1))^2)^0.5;

% Calculate calibration constant
answer      = inputdlg('Distance in mm:');
dist_units  = str2num(answer{1})/1000;

calconst = dist_units ./ dist_pix;

warning on


%--------------------------------------------------------------------------
function isThere = isFileThere(fName,pathh)
fName = [fName '.mat'];
isThere = 0;
a = dir(pathh);
for i = 1:length(a)
    tempName    = a(i).name;
    if length(tempName)==length(fName) && ...
            strcmp(tempName,fName)
        isThere = 1;
    end
end


%--------------------------------------------------------------------------   
function [x,y] = choosePoints(img,link)
%Used for finding coordinate points on a static image 'img'.

warning off

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
hold off

warning on


%--------------------------------------------------------------------------
function [fileName,filePath,fileExt] = giveFileParts(fileName,pathh)

tmp = find(fileName=='.');

if length(tmp)>1
    error('You cannot have periods in the filename');
elseif length(tmp)<1
    error('Your filename needs an extension')
end

fileExt  = fileName(tmp:end);
fileName = fileName(1:tmp-1);
filePath = [pathh filesep fileName fileExt];


%--------------------------------------------------------------------------
function calConstant = zCalibrateM(im)
f       = figure;
set(f,'Color','k');
set(f,'DoubleBuffer','on');

tVal = 0.5;n=0;i=1;x=[];y=[];
while 1 == 1
    imshow(im);hold on
    h = title(['Frame number = ' num2str(i)]);
    set(h,'Color','w');
    plot(x,y,'r+-')
    [xi,yi,but] = ginput(1);
    if isempty(but)
        break
    elseif but==30
        %increase frame
        if i == numFrames,
            i = numFrames;
        else
            i=i+1;
        end
        x = [];y = [];
    elseif but==46
        %increase frames by 10
        if i+10>numFrames
            i = numFrames;
        else 
            i = i+10;
        end
    elseif but==44
        %decrease frames by 10
        if i-10 <1
            i=1;
        else
            i = i-10;
        end
    elseif but==31
        %decrease frame
        if i==1
            i = 1;
        else
            i = i-1;
        end
        x = [];y = [];
    elseif but==102 | but==70
        %Choose first frame
        i = 1;
    elseif but==109 | but==77
        %Choose middle frame
        i = round(numFrames/2);
    elseif but==108 | but==76
        %Choose middle frame
        i = numFrames;
    elseif but==1
        n = n+1;
        x(n) = xi;
        y(n) = yi;
    elseif but==3
        n = n-1;
        x = x(1:n);
        y = y(1:n);
    end
end
close
if ~length(x)==2
    error(['Incorrect number of points!']);
end
dist        = (diff(x).^2 + diff(y).^2).^.5;
disp(' ');
knownDist   = inputdlg('What is the known distance (note: dont give units)?  ');

calConstant = str2num(knownDist{1}) ./ dist;



