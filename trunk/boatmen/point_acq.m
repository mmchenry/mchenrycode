function point_acq(imPath,f_name)
% Acquires landmark data from video for kinematic analysis. Requires that
% two body points arecollected for each frame to define a local coordinate 
% system.
%
% imPath - path to video
% f_name - filename for first frame in video
%
% NOTE: Video must be save as a series of tiff images with filenames that
% end in the frame number.  Check parameter values below for your video.


%% Parameters

% Number of digits at end of image filenames
num_digits = 6; 

% Characters found within each filename that preceed the digits
name_prefix = 'seq';

% Filename suffix
name_suffix = 'tif';


%% Define directories

% Prompt for first frame, if not given
if nargin < 2  
    [f_name,imPath,fIndex] = uigetfile(['*.' name_suffix],...
        'Choose first image in sequence');
    if ~fIndex
        return
    end
end

% Determine number of frames, misc file information
[tPath,tName,tExt,tVers] = fileparts([imPath filesep f_name]);
prefix    = tName(1:end-num_digits);
dirOutput = dir([imPath filesep prefix '*' tExt]);
fileNames = {dirOutput.name}';
numFrames = numel(fileNames);

% Clear variables
clear dirOutput prefix tPath tName tExt tVers


%% Acquire sequence information

% Look for existing seq_info file
a = dir([imPath filesep 'seq_info.mat']);

if isempty(a)
    
    % Prompt user for input
    prompt = {'Individual number','Frame rate (fps)'};
    defaults = {fileNames{1}(7),'30','2'};
    answer = inputdlg(prompt,'Input sequence info',1,defaults);
    
    seq.indiv      = str2num(answer{1});
    seq.frame_rate = str2num(answer{2});
    seq.fileNames  = fileNames;
    seq.numFrames  = numFrames;

    save([imPath filesep 'seq_info'],'seq');
    
    %clear variables
    clear prompt defaults answer
else
    load([imPath filesep 'seq_info']);
end

% Clear variables
clear fileNames numFrames a


%% Look for body coordinate data

% Look for body data
aBody = dir([imPath filesep 'body_data.mat']);

% If no body data . . .
if isempty(aBody)
    
    beep
    disp('You first need to collect body coordinates . . .')
    disp(' ')
    
    % Create structure for body coordinates
    body.headX  = nan(seq.numFrames,1);
    body.headY  = nan(seq.numFrames,1);
    body.tailX  = nan(seq.numFrames,1);
    body.tailY  = nan(seq.numFrames,1);
    
    % Create empty lists for point numbers
    pl(1).point_num = [];
    plNums{1}       = ' ';
    
    % Initialize current frame number
    cFrame = 1;
    %iPoint = 1;
    
else
    % Load body data (structure 'body')
    load([imPath filesep 'body_data.mat'])
    
end


%% Look for appendage data

% Run if body data acquired
if 1 %~isempty(aBody)
    
    % Look for point data
    a = dir([imPath filesep 'appendage_data.mat']);
    
    % Load if file is there
    if ~isempty(a)
        load([imPath filesep 'appendage_data.mat']);
        
        % Make point list
        if isempty(pl)
            plNums = [];
        else
            for i = 1:length(pl)
                plNums{i} = num2str(pl(i).point_num);
            end
        end
    
    % Otherwise, create empty lists
    else
        pl(1).point_num = [];
        plNums{1}       = ' ';
        
    end
    
    
    % Ask what to do this time
    but = questdlg('What do you want to do?','Session question',...
        'Start tracking new point',...
        'Continue an old point',...
        'Delete data on a point',...
        'Start tracking new point');
    
    % Break, if cancelled
    if isempty(but)
        return
    end
    
    switch but
    case 'Start tracking new point'
            
        if isempty(a)
            iPoint = 1;
        else
            iPoint = length(pl)+1;
        end
        
        answer = inputdlg({'What is the point number?'},' ',1,{'1'});
        
        for i = 1:length(plNums)
            if strcmp(plNums{i},answer)
                error('Point number already started')
            end
        end
        
        % Set current frame to one
        cFrame = 1;
        
        % Create fields for new comb plate
        pl(iPoint).point_num = str2num(answer{1});
        pl(iPoint).ptX  = nan(seq.numFrames,1);
        pl(iPoint).ptY  = nan(seq.numFrames,1);
        
        clear answer
            
            
    case 'Continue an old point'
            
        if isempty(a) || isempty(pl)
            error('No plate data exist')
        end
        
        [cFrame,iPoint,pl] = continue_plate(pl,seq);
            
            
    case 'Delete data on a point'
            
        if isempty(a)
            error('No data file exists')
        end
        
        % Prompt for selection
        [s,v] = listdlg('PromptString','Select point number',...
            'SelectionMode','single',...
            'ListString',plNums);
        if isempty(s)
            return
        end
        
        for i=1:length(pl)
            if pl(i).point_num==str2num(plNums{s})
                iPoint = i;
                break
            end
        end
        
        % Confirm delete
        but2 = questdlg(['Are you sure you want to delete point ' ...
            plNums{s} '?'],'Warning','Yes - delete',...
            'No - Cancel','No - Cancel');
        
        % Cancel or execute
        switch but2
        case 'No - Cancel'
            return
                
        case 'Yes - delete'
            if (s==1) && (length(pl)==1)
                pl = [];
            elseif (s==1) && (length(pl)>1)
                pl = pl(2:end);
            elseif s==length(pl)
                pl = pl(1:end-1);
            else
                pl = pl([1:s-1 s+1:end]);
            end
            
            save([imPath filesep 'appendage_data.mat'],'pl')
            return
        end
            
        case 'Cancel' 
            return
    end
    
end



%% Prep for acquisition mode

baseModeColor  = [.8 .4 .4];
angleModeColor = [.7 .7 .9];
mouthModeColor = [.7 .7 0];
statModeColor  = [0 .7 .7];

frameSkip = 5;

% Give instructions
disp(' '); 
disp('=================== CONTROL ====================')
disp('  Left mouse   - pick point'); 
disp('  Right mouse  - remove point'); 
disp('  space bar    - advance 1 frame'); 
disp('  left arrow   - back 1 frame'); 
disp(' ')
disp(['  + - skip ' num2str(frameSkip) ' frames forward']); 
disp(['  - - skip ' num2str(frameSkip) ' frames backward']); 
disp( '  s - change interval for skipping forward and backward');
disp( '  j - jump to frame number');
disp(' ')
disp( '  p - jump to another point');
disp( '  d - toggle displaying all points');
disp( '  c - Toggle contrast display mode'); 
disp( '  z - Zoom mode (return to exit)'); 
disp(' ')
disp( '  t - Point mode');
disp( '  m - Tail mode');
disp( '  o - Head mode');
disp(' ')
disp('Press return or esc when done collecting.')
disp('===============================================')
disp(' ');

% Make figure window
hF = figure;
set(gcf,'DoubleBuffer','on');

% Set initial parameter values
numLimit = 10^10;
progMode = 1;
colorMode = 0;
displayAll = 0;
mouth = 1;

%cFrame = 301


% Display image to get pixel coordinates
warning off
im  = imread([imPath filesep seq.fileNames{cFrame}]); 
hIm = imshow(im);
warning on

% Get x and y limits
xlim_c = xlim;
ylim_c = ylim;


%% Start loop for each new frame
while 1
    
    warning off
    im = imread([imPath filesep seq.fileNames{cFrame}]);
    
    % Convert to grayscale
    if length(size(im))==3
        rgb2gray(im)
    end
    
    % Enhance contrast
    if colorMode
        im = adapthisteq(im,'clipLimit',0.02,...
            'Distribution','rayleigh');
    end
    
    hIm = imshow(im);
    warning on
    
    % Set limits
    xlim(xlim_c);
    ylim(ylim_c);
    
    % Initialize coordinate system
    S = nan(2);
    origin = nan(1,2);
    
    % Change display, depending on mode
    if progMode == 1
        hT = title(['Frame ' num2str(cFrame) ' of ' num2str(seq.numFrames)...
            '  POINT MODE']);
        set(hF,'Color','w')
    elseif progMode == 3
        hT = title(['Frame ' num2str(cFrame) ' of ' num2str(seq.numFrames)...
            '  TAIL MODE']);
        set(hF,'Color',mouthModeColor);
    elseif progMode == 4
        hT = title(['Frame ' num2str(cFrame) ' of ' num2str(seq.numFrames)...
            '  HEAD MODE']);
        set(hF,'Color',statModeColor);
    end
    
    
    %% Start loop for acquiring points
    while 1 
                
        % Plot existing data
        hold on
        h = plotData(pl,cFrame,iPoint,displayAll,body,origin,S);
        hold off
        
        % Prompt for single click
        [x,y,but] = ginput(1);
        
        % If return
        if isempty(but)
            break
            
        % Left click --------------------------------------------------
        elseif but==1 
            
            % Point mode
            if progMode == 1
                pl(iPoint).ptX(cFrame) = x;
                pl(iPoint).ptY(cFrame) = y;
                
            
            % Tail mode
            elseif progMode == 3 
                    body.tailX(cFrame) = x;
                    body.tailY(cFrame) = y;
            
            % Head mode
            elseif progMode == 4 
                    body.headX(cFrame) = x;
                    body.headY(cFrame) = y;        
            end
            
            % Advance frame
            %cFrame = min([cFrame+1 seq.numFrames]);
            break
            
        % Right click --------------------------------------------------  
        elseif but==3 
            
            % Point mode
            if progMode==1
                pl(iPoint).ptX(cFrame) = nan;
                pl(iPoint).ptY(cFrame) = nan;
              
            % Tail mode
            elseif progMode == 3 
                    body.tailX(cFrame) = nan;
                    body.tailY(cFrame) = nan;
                    
            % Head mode        
            elseif progMode == 4 
                    body.headX(cFrame) = nan;
                    body.headY(cFrame) = nan;        
               
            end
            
            break
            
        % If escape    
        elseif but==27 
            x = [];
            y = [];
            break
            
        % If 'z' (zoom)    
        elseif but==122 
            zoom on
            pause;
            xlim_c = xlim;
            ylim_c = ylim;              
            
        % If 'c'
        elseif but==99
            colorMode = abs(colorMode - 1);
            break
            
        % If 'j'
        elseif but==106
            if progMode==2
               beep; disp('You can only look at key frames on angle mode')
            end
            
            answer = inputdlg('Desired frame number',' ',1,{''});
            if ~isempty(answer) && ~isempty(answer{1})
                cFrame = min([seq.numFrames str2num(answer{1})]);
            end
            break
            
         % If 's'
        elseif but==115
            answer = inputdlg('Desired interval',' ',1,{''});
            if ~isempty(answer) && ~isempty(answer{1})
                frameSkip = str2num(answer{1});
            end
            break
            
        % If 'd'
        elseif but==100
            displayAll = abs(displayAll-1);            
        
        % If spacebar or right arrow    
        elseif (but==32) || (but==29) 
            
            mouth = 1;
                
            % If angle mode, advace to next keyframe    
            if progMode==2
                % Locate slices after present that have been modified
                tmp = find(pl(iPoint).baseMod);
                idx = tmp((tmp>cFrame));
                if isempty(idx)
                    beep
                    disp('That is all the keyframes for this plate')
                else
                    cFrame = idx(1);
                end
                
                clear tmp idx
    
            % Advance frame, if in base or tip mode
            else 
                cFrame = min([cFrame+1 seq.numFrames]); 
                
            end
            
            break
            
        % If left arrow
        elseif but== 28
            
            if progMode==2
                
                % Locate slices before present that have been modified
                tmp = find(pl(iPoint).baseMod);
                idx = tmp((tmp<cFrame));
                if isempty(idx)
                    beep
                    disp('No prior keyframes for this plate')
                else
                    cFrame = idx(end);
                end
                
                clear tmp idx
                
            else 
                cFrame = max([cFrame-1 1]);

            end
            
            break
            
        % If -
        elseif but==45
            if (progMode==1) || (progMode==0)
                cFrame = max([cFrame-frameSkip 1]);
                break
            end
            
        % If +   
        elseif but==43 || but==61     
            if (progMode==1) || (progMode==0)
                cFrame = min([cFrame+frameSkip seq.numFrames]);
                break
            end          
        
        % If 'p' (jump to another comb plate)
        elseif but==112
            [cFrame,iPoint,pl] = continue_plate(pl); 
            break
        
        % If 'm' (Tail mode)
        elseif but==109
            
            progMode = 3;
            
            break
        
        % If 'o' (Head mode)
        elseif but==111
            
            progMode = 4;
            
            break  
            
        % If 't' (point mode)    
        elseif but==116
            
            progMode = 1;
            
            break  
   
        end
        
        delete(h);
        
    end
    
    % Save data before changing frames
    save([imPath filesep 'appendage_data.mat'],'pl')
    save([imPath filesep 'body_data.mat'],'body')
    
    % Break out of frame change loop, if return or escape
    if isempty(but) || (but==27)
        break
    end
    
end
close;



function h = plotData(pl,cFrame,iPoint,displayAll,body,origin,S)

% Offset (in ix) for text
offset = 7;

% Define current plate
iPoint_c = iPoint;

hold on

% Include all plates, if displayAll
if displayAll
    iPoint = 1:length(pl);
end

% Plot body points
h(1,1) = plot(body.tailX(cFrame),body.tailY(cFrame),'mo');
h(2,1) = plot(body.headX(cFrame),body.headY(cFrame),'m+');
h(3,1) = plot([body.tailX(cFrame) body.headX(cFrame)],...
            [body.tailY(cFrame) body.headY(cFrame)],'m-');
h(4,1) = plot(origin(1),origin(2),'mo');

% Loop through plates to be displayed
for i = 1:length(iPoint)
    
    % Add fields 'angleX' and 'angleY', if not present in 'pl'
    if ~isfield(pl(iPoint(i)),'angleX') || ...
            isempty(pl(iPoint(i)).angleX)
        pl(iPoint(i)).angleX= nan(length(pl(iPoint(i)).ptX),1);
        pl(iPoint(i)).angleY= nan(length(pl(iPoint(i)).ptX),1);
    end

    % Define points for current plate
    tX = pl(iPoint(i)).ptX(cFrame);
    tY = pl(iPoint(i)).ptY(cFrame);
    
    % Put emphasis on the current plate
    if iPoint_c == iPoint(i)
        tmp = plot(tX,tY,'r+');
    else
        tmp = plot(tX,tY,'w+');
    end
    
    % Add tmp to handle list, clear tmp
    h = [h;tmp]; 
   
    clear tX tY num tmp
end


function [cFrame,iPoint,pl] = continue_plate(pl,seq)

% Make plate list
for i = 1:length(pl)
    plNums{i} = num2str(pl(i).point_num);
end

% Prompt for selection of cureent plate number
[s,v] = listdlg('PromptString','Select plate number',...
    'SelectionMode','single',...
    'ListString',plNums);
if isempty(s)
    return
end

for i=1:length(pl)
    if pl(i).point_num==str2num(plNums{s})
        iPoint = i;
        break
    end
end

% Define curent fraem from end of data
cFrame = find(~isnan(pl(iPoint).ptX),1,'last');
if isempty(cFrame)
    cFrame = 1;
end

clear s v

% Add angleX and angleY, if not present
if ~isfield(pl(iPoint),'angleX') || ...
        isempty(pl(iPoint).angleX)
    pl(iPoint).angleX= nan(seq.numFrames,1);
    pl(iPoint).angleY= nan(seq.numFrames,1);
end


function S = localSystem(P1,P2)
% Defines a transformation vector for a local coordinate system in an
% inertial frame of reference.  Uses P1 as the origin and P2 to find the
% direction of the y-axis.  Coordinates must be (1x2) vectors.

if size(P1,1)~=1 || size(P1,2)~=2 ||...
   size(P2,1)~=1 || size(P2,2)~=2
    error('Coordinates must be (1x2) vectors');
end

yAxis       = (P2-P1)./norm(P2-P1);
xAxis       = [yAxis(2); -yAxis(1)];
S           = [xAxis yAxis'];


function [x,y] = localToGlobal(pts,origin,S)

if size(pts,2)~=2 || size(origin,2)~=2 
    error('Coordinates must be a (nx2) vector');
end

pts         = [inv(S)'*pts']';
pts(:,1)    = pts(:,1)+origin(1);
pts(:,2)    = pts(:,2)+origin(2);
x           = pts(:,1);
y           = pts(:,2);


function [x,y] = globalToLocal(pts,origin,S)

if size(pts,2)~=2 || size(origin,2)~=2 
    error('Coordinates must be a (nx2) vector');
end

pts(:,1)    = pts(:,1)-origin(1);
pts(:,2)    = pts(:,2)-origin(2);
pts         = [S'*pts']';
x           = pts(:,1);
y           = pts(:,2);


