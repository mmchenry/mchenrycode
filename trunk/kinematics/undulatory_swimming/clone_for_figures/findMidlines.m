function mid = findMidlines(vid_dir,bLength,d,roi)
% Finds the body midline for a movie comparised of a series of tif images
%
% Inputs:
% vid_dir - directory containing tifs of video 
% bLength - (1x1) body length of fish (in pixels)
%
% d - structure of head kinematics, including these fields
%   xTip - (nx1) x-coordinates of rostrum
%   yTip - (nx1) y-coordinates of rostrum
%   xCtr - (nx1) x-coordinates of head center
%   yCtr - (nx1) y-coordinates of head center
%   filename - (1xn) string for the filename of tiff image analyzed
%
% roi - (1x2) coordinates for a region of interest used when acquiring head
%             data (optional)
%
% The video frames are assumed to be saved as a tiff series with the final
% digits corresponding to the frame number
% 

% Uses controlfish2_10_c001s0001003876.tif  as the basis for a figure

%% Parameters

visData       = 0;
visSteps      = 1;
%segLength       = (bLength - norm([cntrPnt-rost]))./10;
numEllipsePts = 1000;

edgeThold   = -0.6;


%% File information

files = dir([vid_dir filesep '*.tif']);

% Check match of files and d
if length(d.xTip)~=length(files)
    error('Head data does not match number of tif files')
end

% Find the header at the beginning of file names
name = files(1).name;
iZeros = find(name=='0');
if max(diff(find(name=='0'))) > 1
    firstZero = (iZeros(max(find(diff(iZeros)>1))+1));
else
    firstZero = find(name=='0',1,'first');
end
iNum        = firstZero:length(name);

nameHead    = name(1:firstZero-1);

clear iNum iZeros firstZero name

% Define roi correction, if not given
if nargin < 4
    roi = [0 0];
end


%% Create figure
if visData || visSteps
    hF = figure;
    set(hF,'DoubleBuffer','on')
end

%% Initiate loop through frames

for i = 1:length(files) %397
    
    % Get current filename
    fName = files(i).name;
    
    % Check match with data
    if ~strcmp(fName,d.filename{i})
        error('Video file and data mismatch');
    end
    
    % Get current frame number
    frNum = str2num(files(i).name(length(nameHead)+1:end-4));
    
    % Get image data
    im = imread([vid_dir filesep fName]);
    
    rost            = [d.xTip(i)+roi(1) d.yTip(i)+roi(2)];
    cntrPnt         = [d.xCtr(i)+roi(1) d.yCtr(i)+roi(2)];
    
    % Define points to create initial coordinate system
    P1              = rost;
    P2              = cntrPnt;
    midPnt          = P2;
    bPosition       = norm([P2-P1])/bLength;
    
    % Approximate body width
    bWidth          = norm([rost-cntrPnt])/2;
    lastWidth       = bWidth;
    %maxCent         = -1;
    
    %Lateral dimension of the ellipse
    b     = 2*bWidth(end); 
    
    % Anterio-posterior dimension of ellipse
    a     = b;
    
    % Radial coordinate for ellipse
    theta = linspace(-pi/4,pi/4,numEllipsePts);
    
    % Record data:  
    mid(i).left(1,:)    = P1;
    mid(i).right(1,:)   = P1;
    mid(i).mid(1,:)     = P1;
    mid(i).s(1,1)       = 0;
    mid(i).fName        = fName;
    

%% Step down the body, picking off points 
    
    j = 2;
    
    while 1==1
        
        % 1. Find pixel values along ellipse ------------------------------ 
        
        % Define local coordinate system
        S           = localSystem(P1,P2); 
        
        % Define coordinates for a semi-ellipse        
        circPath(:,1) = a.*cos(theta);
        circPath(:,2) = b.*sin(theta);
        %circPath    = giveHalfEllipse(a,b,numEllipsePts);
        
        % Shift ellipse anteriorly
        circPath    = [circPath(:,1) - 0.5*a circPath(:,2)];
        
        % Rotate to local coordinate system
        circPath    = localToGlobal(circPath,P2,S);
        
        % Interpolate for pixel values along ellipse
        lineVal     = interp2(double(im),circPath(:,1),circPath(:,2));
        
        % Define arclength
        lineDist    = [0 ; cumsum(  (diff(circPath(:,1)).^2 + diff(circPath(:,2)).^2).^.5  )];
        
        % Check if still on frame (i.e. if more than a quarter of pixels
        % are equal to zero)
        if (sum(floor(abs(lineVal))==0)/length(lineVal)) > 0.25 || ...
            (sum(isnan(lineVal))/length(lineVal)) > 0.25
            warning on
            warning(['fell off frame on file ' num2str(i) ': ' fName]);
            break
        end
        
        % Replace nans (from falling off frame) with max values
        if max(isnan(lineVal))==1
            lineVal(isnan(lineVal)) = lineVal(max(find(~isnan(lineVal))));
        end
             
        % 2. Find landmarks along ellipse ---------------------------------
        
        % Normalize pixel values
        lineVal     = (lineVal-max(lineVal))/range(lineVal);
        
        % Low-pass filter the data
        cutFreq = 0.1;
        [B A] = butter(2,cutFreq,'low');
        lineValF = filtfilt(B,A,lineVal);
        
        clear cutFreq B A
        
        % Find indicies for left and right edges
        iLedge = find(lineValF<edgeThold,1,'first');
        iRedge = find(lineValF<edgeThold,1,'last');
        
        % Find index for midpoint
        iMid = round(mean([iLedge iRedge]));
        
        % Define coordinates
        midCoord  = circPath(iMid,:);
        lEdge     = circPath(iLedge,:);
        rEdge     = circPath(iRedge,:);
        bodyWidth = ((lEdge(1)-rEdge(1))^2 + (lEdge(2)-rEdge(2))^2)^0.5;
        
        % Calculate signal noise to determine if still on body
        cutFreq = 0.2;
        [B A] = butter(2,cutFreq,'high');
        noiseLevel = std(filtfilt(B,A,lineVal));
        
        % Normalize noise relative to first iteration
        if j ==2
            noiseLevel0 = noiseLevel;
        end
        noiseLevel = noiseLevel/noiseLevel0;
        
        %Calculate relative body position
        ds = (((mid(i).mid(end,1)-midCoord(1)).^2+...
               (mid(i).mid(end,2)-midCoord(2)).^2).^0.5)./bLength;
        s  = mid(i).s(end,1) + ds;
        
        % Determine whether to break before storing data
        if j>1 &&((s>1.2) || ((s>0.7) && (noiseLevel>2)))
            break
        end
        
        % Visualize results for current iteration
        if 0 %visSteps,  
            warning off
            %subplot(2,1,1)
%             imshow(im);
%             hold on;
%             plot(circPath(:,1),circPath(:,2),'k');
%             plot(midCoord(1),midCoord(2),'r+')
%             plot(lEdge(1),lEdge(2),'g+',rEdge(1),rEdge(2),'y+')
%             title(['bPosition =' num2str(bPosition(end))]); 
%             hold off;
%             subplot(2,1,2);
            plot(lineVal+1,lineDist,'k-')
           % hold on 
%             plot(lineDist,lineValF,'r-')
            plot([lineDist(iLedge) lineDist(iLedge)],[0 -1],'g')
            plot([lineDist(iRedge) lineDist(iRedge)],[0 -1],'y')
            %hold off
            pause(.2)
            warning on
        end      
        
        % Record data:  
        mid(i).left(j,:)    = lEdge;
        mid(i).right(j,:)   = rEdge;
        mid(i).mid(j,:)     = [midCoord];
        mid(i).s(j,1)       = s;
        
        % Set up for next iteration
        P1              = P2;
        P2              = midCoord;
        j               = j+1;       
    
        clear iMid iLedge iRedge midCoord s
    end
    
    if 1 %visData
        warning off
        imshow(im)
        hold on
        plot(mid(i).mid(:,1),mid(i).mid(:,2),'r-o')
        title(['Frame ' num2str(frNum)])
        hold off
        pause(.01)
        warning on
    end
end

%% Close figure
if visData || visSteps
    close(hF)
end


    
function S = localSystem(P1,P2)
%P1 is the origin
    xAxis       = (P2-P1)./norm(P2-P1);
    yAxis       = [xAxis(2);-xAxis(1)];
    S           = [xAxis' yAxis];

function pts = localToGlobal(pts,origin,S)
pts         = [inv(S)'*pts']';
pts(:,1)    = pts(:,1)+origin(1);
pts(:,2)    = pts(:,2)+origin(2);

function pts = globalToLocal(pts,origin,S)
pts(:,1)    = pts(:,1)-origin(1);
pts(:,2)    = pts(:,2)-origin(2);
pts         = [S'*pts']';
