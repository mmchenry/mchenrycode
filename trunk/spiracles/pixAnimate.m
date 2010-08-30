function pixAnimate(pName)
% Animates data collected with pixAcq
%
% pName - path to images and data file
%
% Note: requires that the data file created by pixAcq be in the same
% directory as the image files

%% Parameters

% Name of data file
fName = 'pixel_data.mat';

% Number of times to loop through the data
numLoops = 1;

% Pause duration between screen updates
pauseDur = 0.002;

%% Get path of data file, load data

% if nargin < 1
%     [tmp,pName] = uigetfile(fName,'Select data file');
%     clear tmp
% end

%pName = '/Volumes/Docs/Projects/spiracles/sample sequence/29-6-10 video3';
%pName = 'C:\Users\Guest\Documents\spir\recordings\R006-sp003-r001-29-Jul-2010\pike';
pName = '/Volumes/data_commuter/Projects/Roaches/R005-sp005-r001-29-Jul-2010/pike';


load([pName filesep fName])

%% Make figure

f = figure;
set(f,'DoubleBuffer','on')

subplot(2,1,1)
plot(d.frameNum./d.frameRate,d.pixVal)
xlabel('Time (s)')
ylabel('Mean pixel intensity')
hold on
yVals = ylim;

%% Loop through frames

for j = 1:numLoops
    for i = 1:length(d.frameNum)
        
        % Define current time
        t = d.frameNum(i)./d.frameRate;
        
        % Read current frame
        I = imread([pName filesep d.fileName{i}]);
        
        % Mark current time on graph
        figure(f)
        subplot(2,1,1)
        h1 = plot([t t],yVals,'r-');
        
        % Display current video frame
        subplot(2,1,2);
        warning off
        imshow(I);
        hold on
        warning on
        plot(d.roiX,d.roiY,'r-')
        hold off
        
        % Pause, then delete line
        pause(pauseDur)
        delete(h1)
        
        % Clear variables
        clear t
        
    end
end


