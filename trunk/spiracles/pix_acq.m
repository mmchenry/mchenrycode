function pix_acq(fName,pName,frRate)
% Acquires pixel intensity for a circular region of interest
%
% fName - filename of first file in image sequence
% pName - path for image sequence
% frRate - Frame rate


%% Parameter values
% Header for file names
nameHeader = 'Frame_'; 

% Filename suffix
sfx = 'TIF';

% Number of points for drawing the roi
numRoiPts = 100;

% Number of decimal places for frame numbers filename
numDigits = 4;

% Show the acquisition
visSteps = 0;

% Window size for x-correlation analysis
winSize = 256;

% Interval for advancing frames
frameSkip = 1;


%% Initialize for recording

% Prompt to select image sequence
if nargin < 2
    [fName,pName] = uigetfile(['*.' sfx ],'Select first image in sequence');
    if pName==0
        return
    end
end

% Prompt for frame rate
if nargin < 3
    answer = inputdlg('Frame rate (fps):',' ',1,{'15'});
    if isempty(answer)
        return
    end
    frRate = str2num(answer{1});
end

% Check name header
if ~strcmp(nameHeader,fName(1:length(nameHeader)))
    error('Name header incorrect')
end

% Check number of frames in directory
a = dir([pName filesep '*.' sfx]);
numFrames = length(a);

% Store parameter values
d.frameRate = frRate;
d.numFrames = numFrames;
d.path      = pName;
d.nameHead  = nameHeader;
d.numDigits = numDigits;

clear frameRate


%% Interactively define roi

% Show first frame in figure window
f = figure;
I = imread([pName filesep a(1).name]);
imshow(I)
hold on

% Select center point
title('Choose center point')
[xStart,yStart] = ginput(1);

% Store center for first frame
d.xCntr = xStart;
d.yCntr = yStart;

% Plot center
plot(d.xCntr,d.yCntr,'r+')

% Define unit circle
theta = linspace(0,2*pi,numRoiPts);
d.xC = cos(theta);
d.yC = sin(theta);

% Select circular roi
title('Choose radius')
[xT,yT] = ginput(1);
d.r = ((d.xCntr-xT)^2 + (d.yCntr-yT)^2)^0.5;
xROI = d.r.*d.xC + d.xCntr;
yROI = d.r.*d.yC + d.yCntr;

% Binary at roi
roiI0 = roipoly(I,xROI,yROI);

% Show selected roi
plot(xROI,yROI,'r-')
pause(0.5)
close


clear theta


%% Load initial image

% % Define unit circle
% theta = linspace(0,2*pi,numRoiPts);
% xCirc = cos(theta);
% yCirc = sin(theta);

% Current frame number
frNum = ['0001'];
frName = [nameHeader frNum((end-numDigits+1):end) '.' sfx];
I0 = imread([pName filesep frName]);

% Interrogation window
% winI0 = I0(floor(d.yCntr(1)-winSize/2):floor(d.yCntr(1)+winSize/2),...
%              floor(d.xCntr(1)-winSize/2):floor(d.xCntr(1)+winSize/2));

%winSize = min([size(I0,1) size(I0,2)]);

winI0 = I0(ceil(size(I0,1)/2-winSize/2+1):floor(size(I0,1)/2+winSize/2-1),...
           ceil(size(I0,2)/2-winSize/2+1):floor(size(I0,2)/2+winSize/2-1));
         
clear frName frNum 


%% Track frame displacement with cross-correlation

if visSteps
    hf = figure;
    set(hf,'DoubleBuffer','on')
    pause(.1)
end

% Maximum displacement allowed before resetting the 
maxDisp = 10;

% Check for file
a = dir([pName filesep 'roi_data.mat']);

if isempty(a)
    
    d.frames = 1:frameSkip:numFrames;
    
    % Loop through frames
    for i = 2:length(d.frames)
        
        % Current frame number
        frNum = ['000' num2str(d.frames(i))];
        frName = [nameHeader frNum((end-numDigits+1):end) '.' sfx];
        
        % Read image
        I = imread([pName filesep frName]);
        
        % Define roi in image
        winI = I(ceil(size(I0,1)/2-winSize/2+1):floor(size(I0,1)/2+winSize/2-1),...
            ceil(size(I0,2)/2-winSize/2+1):floor(size(I0,2)/2+winSize/2-1));
        
        % Find offset from image window pair
        %[xOffset,yOffset,totOffset] = findOffset(winI,winI0);
        
        % Correlation coefficient for interrogation windows
        cc = normxcorr2(winI,winI0);
        
        % Find index of max cc value
        [max_cc, imax] = max(cc(:));
        
        % Find coordinates from max index
        [ypeak, xpeak] = ind2sub(size(cc),imax(1));
        
        % Determine offset from coordinates
        yOffset = ypeak-size(winI0,1);
        xOffset = xpeak-size(winI0,2);
        
        % Find total displacement
        totOffset = sqrt(xOffset^2 + yOffset^2);
        
        % Update center point from offset
        d.xCntr(i) = xStart - xOffset;
        d.yCntr(i) = yStart - yOffset;
        
        % Define roi coordinates
        xROI = d.r.*d.xC + d.xCntr(i);
        yROI = d.r.*d.yC + d.yCntr(i);
        
        % Define roi binary
        roiI = roipoly(I,xROI,yROI);
        
        % Update status
        disp(['Done frame ' num2str(i) ' of ' num2str(numFrames)])      
        
        if (totOffset > maxDisp) 
   
            % Reset reference image and center point
            disp('Resetting reference . . .')
            xStart = d.xCntr(i);
            yStart = d.yCntr(i);
            I0 = I;
            winI0 = I0(ceil(size(I0,1)/2-winSize/2+1):floor(size(I0,1)/2+winSize/2-1),...
                       ceil(size(I0,2)/2-winSize/2+1):floor(size(I0,2)/2+winSize/2-1));
        end
        
        if visSteps
            % Coordinates for cross-correlation data
            xVals = [-(size(winI0,2)-1):1:size(winI0,2)-1];
            yVals = [-(size(winI0,1)-1):1:size(winI0,1)-1];
            
            warning off
            
            subplot(3,3,1)
            imshow(winI0)
            title('base image')
            
            subplot(3,3,4)
            imshow(winI)
            title('current image')
            
            subplot(3,3,7)
            h = pcolor(xVals,yVals,cc);
            
            axis square
            set(h,'EdgeColor','none');
            hold on
            title(['Total displacement = ' num2str(totOffset)])
            
            plot([0 0],ylim,'k',xlim,[0 0],'k',xOffset,yOffset,'ko');
            hold off
            colormap jet
            
            subplot(3,3,[2 3 5 6 8 9])
            imshow(I)
            hold on
            plot(xROI,yROI,'r-')
            hold off
            
            pause(.001)
            warning on
        end
        
        % Clear variables
        clear I frName frNum winI xOffset yOffset cc max_cc imax totOffset roiI
        clear xROI yROI
    end
    
    clear a maxDisp
    
    save([pName filesep 'roi_data'],'d')
    
else
    
    disp(' ')
    disp('Loading previous roi data file . . .')
    load([pName filesep 'roi_data'])
    
end


%% Postprocessing of roi kinematics

cut_freq = 0.05;

d.xCntr_f = butter_filt(d.xCntr,1,cut_freq,'low'); 
d.yCntr_f = butter_filt(d.yCntr,1,cut_freq,'low'); 

if 0
    figure
    subplot(2,1,1)
    plot(d.frames,d.xCntr,'.',d.frames,d.xCntr_f,'-')
    ylabel('x')
    subplot(2,1,2)
    plot(d.frames,d.xCntr,'.',d.frames,d.xCntr_f,'-')
    ylabel('y')
end


%% Evaluate pixel values


% Loop through frames
for i = 1:length(d.frames)

    % Current frame number
    frNum = ['000' num2str(d.frames(i))];
    frName = [nameHeader frNum((end-numDigits+1):end) '.' sfx];
    
    % Read image
    I = imread([pName filesep frName]);
    
    % Define roi coordinates
    d.xROI{i} = d.r.*d.xC + d.xCntr(i);
    d.yROI{i} = d.r.*d.yC + d.yCntr(i);
    
    % Define roi binary
    roiI = roipoly(I,d.xROI{i},d.yROI{i});
    
    % Visualize roi analyzed
    if visSteps      
        
        %tmp(roiI(:))=I(roiI(:));     
        bI = uint8(255.*ones(size(I,1),size(I,2)));
        bI(roiI(:))=I(roiI(:));
        
        warning off
        subplot(1,2,1)
        imshow(bI)
        axis([min(d.xROI{i}) max(d.xROI{i}) ...
            min(d.yROI{i}) max(d.yROI{i})])
        
        title(['Corrected: Frame ' num2str(d.frames(i)) ' of ' num2str(d.frames(end))]);
            
        bI = uint8(255.*ones(size(I,1),size(I,2)));
        bI(roiI0(:))=I(roiI0(:));
        
        subplot(1,2,2)
        imshow(bI)
        axis([min(d.xROI{1}) max(d.xROI{1}) ...
            min(d.yROI{1}) max(d.yROI{1})])
        title('Uncorrected')
        warning on
        
        pause(.1)
        
    else
        disp(['Measuring pixel intensity: done frame ' ...
              num2str(d.frames(i)) ' of ' num2str(d.frames(end))])
    end
    
    % Calc & store mean pixel value
    d.pixVal(i)    = mean(double(I(roiI(:))));
    d.fileName{i}  = frName;
    d.frameNum(i)  = str2num(frNum);
    
    % Clear variables
    clear I frName frNum winI xOffset yOffset
end

% if ~visSteps
%     figure
% end

% Plot data
% plot(d.frameNum./d.frameRate,d.pixVal)
% xlabel('Time (s)')
% ylabel('Mean pixel intensity')

% Save data
%[sFileName,sFilePath]= uiputfile('pixdata.mat','Save data');
%save([sFilePath filesep sFileName],'d');

save([pName filesep 'pixel_data'],'d')



function [xOffset,yOffset,totOffset] = findOffset(winI,winI0)
% Find offset from image window pair

% Correlation coefficient for interrogation windows
cc = normxcorr2(winI,winI0);

% Find index of max cc value
[max_cc, imax] = max(abs(cc(:)));

% Find coordinates from max index
[ypeak, xpeak] = ind2sub(size(cc),imax(1));

% Determine offset from coordinates
yOffset = ypeak-size(winI0,1);
xOffset = xpeak-size(winI0,2);

% Find total displacement
totOffset = sqrt(xOffset^2 + yOffset^2);

if 1
    % Coordinates for cross-correlation data
    xVals = [-(size(winI0,2)-1):1:size(winI0,2)-1];
    yVals = [-(size(winI0,1)-1):1:size(winI0,1)-1];
    
    subplot(3,1,1)
    imshow(winI0)
    title('base image')
    
    subplot(3,1,2)
    imshow(winI)
    title('current image')
    
    subplot(3,1,3)
    h = pcolor(xVals,yVals,cc);
    
    axis square
    set(h,'EdgeColor','none');
    hold on
    title('cross-correlation')
    plot([0 0],ylim,'k',xlim,[0 0],'k',xOffset,yOffset,'ko');
    hold off
    colormap jet
end

return


 
    %     % Define roi mask
%     xCir = r.*xC + xCntr + xDisp;
%     yCir = r.*yC + yCntr + yDisp;
%     roiI = roipoly(I,xCir,yCir);
%     I(roiI(:))=I(roiI(:));
%     bI(roiI(:))=I(roiI(:));
    
    %roiI = [roiI(:,1)+xDisp roiI(:,1)+yDisp];
    
    % Visualize cross-correlation approach
    if 0 %visSteps
        
        % Correct for yOffset
        if yOffset<0
            I_off = I(abs(yOffset)+1:end,:);
        elseif yOffset==0
            I_off = I;
        else
            I_off = I(1:end-abs(yOffset)-1,:);
        end
        
        % Correct for xOffset
        if xOffset<0
            I_off = I_off(:,abs(xOffset)+1:end);
        else
            I_off = I_off(:,1:end-abs(xOffset)-1);
        end
        % Interrogation window
    winIoff = I_off(floor(d.yCntr(i)-winSize/2):floor(d.yCntr(i)+winSize/2),...
                    floor(d.xCntr(i)-winSize/2):floor(d.xCntr(i)+winSize/2));
        
        
        
        bI = uint8(mean(I(:)).*zeros(size(I,1),size(I,2)));
        
        
        %figure;
        warning off
 
        subplot(2,2,1)
        imshow(winI0)
        title('base image')

        subplot(2,2,3)
        imshow(winI)
        title('current image')
        
  
        
        subplot(2,2,2)
        imshow(winIoff)
        title('corrected current image')  
        
        colormap jet
    end

    
    
function data_filtered = butter_filt(data,sample_rate,cut_freq,type) 
% High-pass or low-pass butterworth filter

% All frequency values are in Hz.

% Nyquist freq.
Nqst = sample_rate/2;   

% Calculate stopband frequency
if strcmp(type,'high')
    stop_freq = max([(cut_freq - Nqst/10) .01]);  

elseif strcmp(type,'low')
    stop_freq = min([(cut_freq + Nqst/10) (Nqst-.01)]); 
 
end

% Stopband Attenuation (dB)
Astop = 30;

% Passband Ripple (dB)
Apass = 1;   

% Normalise the cutoff freq. to the Nyquist freq
Wp = cut_freq/Nqst;

% Normalise the stoppass freq. to the Nyquist freq
Ws = stop_freq/Nqst;

% Check cutoff
if (Wp > 1) || (Ws > 1)
    error('Cutoff or bandpass frequencies cannot exceed the Nyquist freq')
end

% Calculate the order from the parameters using BUTTORD.
[N,Fc] = buttord(Wp, Ws, Apass, Astop);    
    
% Calculate the zpk values using the BUTTER function.
[B A] = butter(N, Fc, type);

% Plot frequency reponse
%freqz(B,A,512,sample_rate); 

% Filter the data
data_filtered   = filtfilt(B,A,data); 

