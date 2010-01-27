function acquire_linkage(calConst,imPath)
% Interactively collects coordinate points from a stomatopod appendage
%
% calConst assumes units of m/pixel.
% rootPath is the root directory for the image directories used in the
% analysis


%% Prompt for image file, if not given

if nargin < 2
    [fileName,pathh] = uigetfile('*.*','Choose an image file');
end

[imName,imPath,imExt] = giveFileParts(fileName,pathh);

clear fileName


%% Determine the calibration constant, if not given

if nargin < 1
    numReplicates       = 1;
    [cfileName,cpathh]    = uigetfile('*.*','Choose a calibration image file');
    [cfileName,cfilePath,cfileExt] = giveFileParts(cfileName,cpathh);
    imgCal              = imread([cpathh filesep cfileName],cfileExt(2:end));
    calData             = runCalibrations(imgCal,numReplicates);
    calConst            = calData.meanCalConst 
end

clear cpathh calData imgCal cfilePath cfileExt imgCal numReplicates


%% Interactively acquire coordinates

img = imread(imPath);

% Text to place next to each coordinate
txt{1} = 'A';
txt{2} = 'B';
txt{3} = 'C';
txt{4} = 'D';
txt{5} = 'E';
txt{6} = 'F';
txt{7} = 'G';
txt{8} = 'H';

% Title text to tell user what to collect, based on Fig. 2 of 
% Patek et al (2007), with addition of E, F, G, & H
ttl{1} = 'Choose A: merus, meryl-V joint';
ttl{2} = 'Choose B: meryl-V, carpus joint';
ttl{3} = 'Choose C: carpus, merus joint';
ttl{4} = 'Choose D: base of saddle';
ttl{5} = 'Choose E: top of dactyl suture';
ttl{6} = 'Choose F: Point of contact on dactyl';
ttl{7} = 'Choose G: Distal end of dactyl';
ttl{8} = 'Choose H: proximal end of merus';
ttl{9} = 'Press return to finish';

% Prep figure
figure;
imshow(img);
title(ttl{1});
set(gcf,'DoubleBuffer','on');

% Give instructions
disp(' '); disp(' ');
disp('Left mouse button picks points.');disp(' ');
disp('Right mouse button removes last point.');disp(' ');
disp('Press return to stop.')

% Loop through point collection
n = 0; but = 1;

while 1 == 1
    [xi,yi,but] = ginput(1);
    
    % Return
    if n==8 && isempty(but)
        break
    % Return & data complete    
    elseif isempty(but)
        error('Not done collecting -- data not saved');
        return

    % Left button
    elseif but==1
        n = min([n+1 8]);
        x(n) = xi;
        y(n) = yi;
        
    % Right button
    elseif but==3
        if n-1 == 0
            n = 0;
            x = [];
            y = [];
        else
            n = n-1;
            x = x(1:n);
            y = y(1:n);
        end
        
    end
    
    % Update image
    imshow(img);
    hold on
    plot(x,y,'ro')
    title(ttl{n+1});
    hold off
    
    if ~(n==0)
        hold on
        for i=1:n
            tt = text(x(i)+15,y(i)+15,txt{i});
            set(tt,'Color','r')
        end
        hold off 
    end
    
end

close 

% Store data in 'coords'
coords.txt      = txt;
coords.info     = ttl;
coords.x_pix    = x;
coords.y_pix    = y;
coords.x        = x .*calConst;
coords.y        = y .*calConst;
coords.calConst = calConst;

% Save in image directory
save([pathh filesep imName '_linkdata'],'coords');




%% FUNCTIONS

%--------------------------------------------------------------------------
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
            i == numFrames;
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

