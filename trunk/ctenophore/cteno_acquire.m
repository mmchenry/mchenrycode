function cteno_acquire(imPath,f_name)


%% Parameters

num_digits = 6;     % Digits at end of image filenames
name_prefix = 'cteno';


%% Define directories

% Prompt for first frame, if not given
if nargin < 2  
    [f_name,imPath,fIndex] = uigetfile({'*.tif','*.jpg'},...
        'Choose first image in sequence');
    if ~fIndex
        return
    end
end

% Determine number of frames, misc file information
[tPath,tName,tExt,tVers] = fileparts([imPath filesep f_name]);
prefix = tName(1:end-num_digits);
dirOutput = dir([imPath filesep prefix '*' tExt]);
fileNames = {dirOutput.name}';
numFrames = numel(fileNames);

clear dirOutput prefix tPath tName tExt tVers


%% Acquire sequence information

% Look for existing seq_info file
a = dir([imPath filesep 'seq_info.mat']);

if isempty(a)
    % Prompt user for input
    prompt = {'Individual number','Frame rate (fps)'};
    defaults = {fileNames{1}(7),'30'};
    answer = inputdlg(prompt,'Input sequence info',1,defaults);
    
    seq.indiv = num2str(answer{1});
    seq.frame_rate = num2str(answer{2});
    seq.fileNames = fileNames;
    seq.numFrames = numFrames;
    
    save([imPath filesep 'seq_info'],'seq');
    
    %clear variables
    clear prompt defaults answer
else
    load([imPath filesep 'seq_info']);
end

% Clear variables
clear fileNames numFrames


%% Prompt for what to do in this session

% Look for body data
aBody = dir([imPath filesep 'body_data.mat']);

% Load or create structure, if no file
if isempty(aBody)

    beep
    disp('You first need to collect body coordinates . . .')
    disp(' ')
    
    body.statX  = nan(seq.numFrames,1);
    body.statY  = nan(seq.numFrames,1);
    body.mouthX = nan(seq.numFrames,1);
    body.mouthY = nan(seq.numFrames,1);
    
else   
    % Load body data (structure 'body')
    load([imPath filesep 'body_data.mat'])
end

% Look for plate data
a = dir([imPath filesep 'plate_data.mat']);

% Load or create structure, if no file
if ~isempty(a)
    
    load([imPath filesep 'plate_data.mat']);
    
    % Make plate list
    if isempty(pl)
        plNums = [];
    else
        for i = 1:length(pl)
            plNums{i} = num2str(pl(i).plate_num);
        end
    end
    
else
    
    % Otherwise, cerate empty lists
    pl(1).plate_num = [];
    plNums{1}       = ' ';
    
end

% Ask what to do this time
but = questdlg('What do you want to do?','Session question',...
    'Start tracking new comb plate',...
    'Continue an old comb plate',...
    'Delete data on a comb plate',...
    'Start tracking new comb plate');

% Break, if cancelled
if isempty(but)
    return
end

switch but
    case 'Start tracking new comb plate'
        
        if isempty(a)
            iPlate = 1;
        else
            iPlate = length(pl)+1;
        end
        
        answer = inputdlg({'What is the plate number?'},' ',1,{'1'});
        
        for i = 1:length(plNums)
            if strcmp(plNums{i},answer)
                error('Plate number already started')
            end
        end
        
        % Set current frame to one
        cFrame = 1;
        
        % Create fields for new comb plate
        pl(iPlate).plate_num = str2num(answer{1});
        pl(iPlate).baseX = nan(seq.numFrames,1);
        pl(iPlate).baseY = nan(seq.numFrames,1);
        pl(iPlate).baseMod = zeros(seq.numFrames,1);
        pl(iPlate).tipX  = nan(seq.numFrames,1);
        pl(iPlate).tipY  = nan(seq.numFrames,1);
        pl(iPlate).angleX= nan(seq.numFrames,1);
        pl(iPlate).angleY= nan(seq.numFrames,1);
        
        clear answer
        
        
    case 'Continue an old comb plate'
        
        if isempty(a) || isempty(pl)
            error('No plate data exist')
        end
        
        [cFrame,iPlate,pl] = continue_plate(pl);
        
        
    case 'Delete data on a comb plate'
        
        if isempty(a)
            error('No data file exists')
        end
        
        
        % Prompt for selection
        [s,v] = listdlg('PromptString','Select plate number',...
            'SelectionMode','single',...
            'ListString',plNums);
        if isempty(s)
            return
        end
        
        for i=1:length(pl)
            if pl(i).plate_num==str2num(plNums{s})
                iPlate = i;
                break
            end
        end
        
        % Confirm delete
        but2 = questdlg(['Are you sure you want to delete plate ' ...
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
                
                save([imPath filesep 'plate_data.mat'],'pl')
                return
        end
        
    case 'Cancel'
        
        return
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
disp( '  p - jump to another comb plate');
disp( '  d - toggle displaying all comb plates');
disp( '  c - Toggle color display mode'); 
disp( '  z - Zoom mode (return to exit)'); 
disp(' ')
disp( '  b - Base mode');
disp( '  t - Tip mode');
disp( '  a - Angle mode');
disp( '  m - Mouth mode');
disp( '  o - Statocyst mode');
disp(' ')
disp( '  r - run body coord analysis');
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
    
    %im = histeq(rgb2gray(im));
    if ~colorMode
        im = adapthisteq(rgb2gray(im),'clipLimit',0.02,...
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
            '  TIP MODE']);
        set(hF,'Color','w')
    elseif progMode == 0
        hT = title(['Frame ' num2str(cFrame) ' of ' num2str(seq.numFrames)...
            '  BASE MODE']);
        set(hF,'Color',baseModeColor);
    elseif progMode == 2
        hT = title(['Frame ' num2str(cFrame) ' of ' num2str(seq.numFrames)...
            '  ANGLE MODE']);
        set(hF,'Color',angleModeColor);
    elseif progMode == 3
        hT = title(['Frame ' num2str(cFrame) ' of ' num2str(seq.numFrames)...
            '  MOUTH MODE']);
        set(hF,'Color',mouthModeColor);
    elseif progMode == 4
        hT = title(['Frame ' num2str(cFrame) ' of ' num2str(seq.numFrames)...
            '  STATOCYST MODE']);
        set(hF,'Color',statModeColor);
    end
    
    
    %% Start loop for acquiring points
    while 1 
                
        % Plot existing data
        hold on
        h = plotData(pl,cFrame,iPlate,displayAll,body,origin,S);
        hold off
        
        % Prompt for single click
        [x,y,but] = ginput(1);
        
        % If return
        if isempty(but)
            break
            
        % Left click --------------------------------------------------
        elseif but==1 
            
            % Tip mode
            if progMode == 1
                pl(iPlate).tipX(cFrame) = x;
                pl(iPlate).tipY(cFrame) = y;
                
                % Base mode
            elseif progMode == 0
                
                % Check that mouth and statocyst points are defined
                if isnan(body.mouthX(cFrame))
                    beep
                    disp(' ')
                    disp('You need to first define the mouth coordinate for this frame')
                    progMode = 3;
                    break
                    
                elseif isnan(body.statX(cFrame))
                    beep
                    disp(' ')
                    disp('You need to define the statocyst coordinate for this frame')
                    progMode = 4;
                    break
                    
                % Otherwise, define coordinate system
                else
                    origin = [mean([body.mouthX(cFrame) body.statX(cFrame)]) ...
                        mean([body.mouthY(cFrame) body.statY(cFrame)])];
                    S = localSystem(origin,[body.mouthX(cFrame) body.mouthY(cFrame)]);
                end
            
                % Set that the current base value has been modified
                pl(iPlate).baseMod(cFrame) = 1;
                
                % Define frames that have body coordinates
                idx_b = ~isnan(body.mouthX) & ~isnan(body.statX);
                
                % Among those having body coordinate, locate keyframes 
                % after present 
                tmp = find(pl(iPlate).baseMod & idx_b);
                idx = tmp((tmp>cFrame));
                
                % Define interval over which to define the base coords
                % (i.e. up to next keyframe)
                if isempty(idx)
                    iFill = min([cFrame seq.numFrames]):find(idx_b,1,'last');
                else
                    iFill = min([cFrame seq.numFrames]):idx(1);
                end
                
                % Define base coordinate in body system
                [xB,yB] = globalToLocal([x y],origin,S);
                
                % Loop through each frame of the interval
                for i = 1:length(iFill)
                    
                    % Check for gap in body data
                    if isnan(body.statX(iFill(i))) || ...
                            isnan(body.mouthX(iFill(i)))
                        break
                    end
                    
                    % Define origin and transform
                    or_c = [mean([body.mouthX(iFill(i)) body.statX(iFill(i))]) ...
                            mean([body.mouthY(iFill(i)) body.statY(iFill(i))])];
                    S_c  = localSystem(origin,[body.mouthX(iFill(i)) ...
                                               body.mouthY(iFill(i))]);
                                           
                    % Transform body system point into global coordinates
                    [x_c,y_c] = localToGlobal([xB yB],or_c,S_c);
                    
                    % Store coordinates
                    pl(iPlate).baseX(iFill(i)) = x_c;
                    pl(iPlate).baseY(iFill(i)) = y_c;
                    
                    clear x_c y_c or_c S_c
                end
                
                clear xB yB iFill idx_b tmp
                
                % Update base coords
                %pl(iPlate).baseX(iFill) = x.*ones(length(iFill),1);
                %pl(iPlate).baseY(iFill) = y.*ones(length(iFill),1);
            
            % Angle mode
            elseif progMode == 2
                % Check that mouth and statocyst points are defined
                if isnan(body.mouthX(cFrame))
                    beep
                    disp(' ')
                    disp('You need to first define the mouth coordinate for this frame')
                    progMode = 3;
                    break
                    
                elseif isnan(body.statX(cFrame))
                    beep
                    disp(' ')
                    disp('You need to define the statocyst coordinate for this frame')
                    progMode = 4;
                    break
                    
                % Otherwise, define coordinate system
                else
                    origin = [mean([body.mouthX(cFrame) body.statX(cFrame)]) ...
                        mean([body.mouthY(cFrame) body.statY(cFrame)])];
                    S = localSystem(origin,[body.mouthX(cFrame) body.mouthY(cFrame)]);
                end
                
                % Define frames that have body coordinates
                idx_b = ~isnan(body.mouthX) & ~isnan(body.statX);
                
                % Among those having body coordinate, locate keyframes 
                % after present 
                tmp = find(pl(iPlate).baseMod & idx_b);
                idx = tmp((tmp>cFrame));
                
                % Define interval over which to define the base coords
                % (i.e. up to next keyframe)
                if isempty(idx)
                    iFill = min([cFrame seq.numFrames]):find(idx_b,1,'last');
                else
                    iFill = min([cFrame seq.numFrames]):idx(1);
                end
                
                % Define base coordinate in body system
                [xB,yB] = globalToLocal([x y],origin,S);
                
                % Loop through each frame of the interval
                for i = 1:length(iFill)
                    
                    % Check for gap in body data
                    if isnan(body.statX(iFill(i))) || ...
                            isnan(body.mouthX(iFill(i)))
                        break
                    end
                    
                    % Define origin and transform
                    or_c = [mean([body.mouthX(iFill(i)) body.statX(iFill(i))]) ...
                            mean([body.mouthY(iFill(i)) body.statY(iFill(i))])];
                    S_c  = localSystem(origin,[body.mouthX(iFill(i)) ...
                                               body.mouthY(iFill(i))]);
                                           
                    % Transform body system point into global coordinates
                    [x_c,y_c] = localToGlobal([xB yB],or_c,S_c);
                    
                    % Store coordinates
                    pl(iPlate).angleX(iFill(i)) = x_c;
                    pl(iPlate).angleY(iFill(i)) = y_c;
                    
                    clear x_c y_c or_c S_c
                end
                
                clear xB yB iFill idx_b tmp
                
                
%                 % Locate slices after present that have been modified
%                 tmp = find(pl(iPlate).baseMod);
%                 idx = tmp((tmp>cFrame));
%                 
%                 % Define interval over which to define the base coords
%                 if isempty(idx)
%                     iFill = min([cFrame seq.numFrames]):seq.numFrames;
%                 else
%                     iFill = min([cFrame seq.numFrames]):idx(1);
%                 end
%                 
%                 % Update angle coords
%                 pl(iPlate).angleX(iFill) = x.*ones(length(iFill),1);
%                 pl(iPlate).angleY(iFill) = y.*ones(length(iFill),1);
            
            % Mouth mode
            elseif progMode == 3 
                    body.mouthX(cFrame) = x;
                    body.mouthY(cFrame) = y;
            
            % Statocyst mode
            elseif progMode == 4 
                    body.statX(cFrame) = x;
                    body.statY(cFrame) = y;        
            end
            
            % Advance frame
            %cFrame = min([cFrame+1 seq.numFrames]);
            break
            
        % Right click --------------------------------------------------  
        elseif but==3 
            
            % Tip mode
            if progMode==1
                pl(iPlate).tipX(cFrame) = nan;
                pl(iPlate).tipY(cFrame) = nan;
            
            % Base mode
            elseif progMode==0
                
                % Locate slices after present that have been modified
                tmp = find(pl(iPlate).baseMod);
                idx = tmp((tmp>cFrame));
                
                % Define interval over which to define the base coords
                if pl(iPlate).baseMod(cFrame)
                    
                    % Set that the current base value has been unmodified
                    pl(iPlate).baseMod(cFrame) = 0;
                    
                    if isempty(idx)
                        iFill = min([cFrame seq.numFrames]):seq.numFrames;
                    else
                        iFill = min([cFrame seq.numFrames]):idx(1);
                    end
                    
                else
                    iFill = cFrame;
                end
                
                % Update base coords
                pl(iPlate).baseX(iFill) = nan(length(iFill),1);
                pl(iPlate).baseY(iFill) = nan(length(iFill),1);
                
                % Update angle coords
                pl(iPlate).angleX(iFill) = nan(length(iFill),1);
                pl(iPlate).angleY(iFill) = nan(length(iFill),1);
            
            % Angle mode
            elseif progMode==2
                
                % Locate slices after present that have been modified
                tmp = find(pl(iPlate).baseMod);
                idx = tmp((tmp>cFrame));
                
                % Define interval over which to define the base coords
                if isempty(idx)
                    iFill = min([cFrame seq.numFrames]):seq.numFrames;
                else
                    iFill = min([cFrame seq.numFrames]):idx(1);
                end
                
                % Update angle coords
                pl(iPlate).angleX(iFill) = nan(length(iFill),1);
                pl(iPlate).angleY(iFill) = nan(length(iFill),1);
                
            elseif progMode == 3 % Mouth mode
                    body.mouthX(cFrame) = nan;
                    body.mouthY(cFrame) = nan;
                    
            elseif progMode == 4 % Statocyst mode
                    body.statX(cFrame) = nan;
                    body.statY(cFrame) = nan;        
               
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
            
        % If 'b'
        elseif (but==98) 
            
           progMode = 0;
           break
            
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
                tmp = find(pl(iPlate).baseMod);
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
                tmp = find(pl(iPlate).baseMod);
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
            
        % If 'a' (enter angle mode)
        elseif but==97
            
            if progMode ~=2
                answer = questdlg(['Enter in this mode only after ' ...
                    'collecting all base and tip data'],'Warning',...
                    'Continue','Cancel','Continue');
                if strcmp(answer,'Continue')
                    
                    % Set progMode
                    progMode = 2;
                    
                    % Jump to first keyframe
                    tmp = find(pl(iPlate).baseMod,1,'first');
                    if isempty(tmp)
                        beep; disp('No keyframes defined')
                        progMode = 0;
                    else
                        cFrame = tmp;
                    end
                    
                    clear tmp
                end
            else
                progMode = 1;
            end
            
            break
        
        % If 'p' (jump to another comb plate)
        elseif but==112
            [cFrame,iPlate,pl] = continue_plate(pl); 
            break
        
        % If 'm' (mouth landmark mode)
        elseif but==109
            
            progMode = 3;
            
            break
        
        % If 'o' (statocyst mode)
        elseif but==111
            
            progMode = 4;
            
            break  
            
        % If 't' (tip mode)    
        elseif but==116
            
            progMode = 1;
            
            break  
        
        % If 'r' (run body coordinate analysis)
        elseif but==114
            answer = questdlg(['Enter in this mode only after ' ...
                'collecting all body coordinates'],'Warning',...
                'Continue','Cancel','Continue');
            
            if strcmp(answer,'Continue')
                
                % Get current title, update title
                disp(' ')
                disp('Running body coordinate analysis . . .')
                
                %disp('This feature has not been implemented yet');
                body = bodyAnalysis(body);
                
            end
            
            disp(' . . . done');
            disp(' ');
            
            break
            
        end
        
        delete(h);
        
    end
    
    % Save data before changing frames
    save([imPath filesep 'plate_data.mat'],'pl')
    save([imPath filesep 'body_data.mat'],'body')
    
    % Break out of frame change loop, if return or escape
    if isempty(but) || (but==27)
        break
    end
    
end
close;

function body = bodyAnalysis(body)

% Extract points
if isfield(body,'raw')
    sX = body.raw.statX;
    sY = body.raw.statY;
    mX = body.raw.mouthX;
    mY = body.raw.mouthY;
    
else
    sX = body.statX;
    sY = body.statY;
    mX = body.mouthX;
    mY = body.mouthY;
    
    % Store data away in raw
    body.raw.statX  = sX;
    body.raw.statY  = sY;
    body.raw.mouthX = mX;
    body.raw.mouthY = mY;
end

% Define index
idx = 1:length(sX);

% Interpolate to remove nans between points
warning off

iNan = isnan(sX);
sX(iNan) = interp1(idx(~iNan),sX(~iNan),idx(iNan));

iNan = isnan(sY);
sY(iNan) = interp1(idx(~iNan),sY(~iNan),idx(iNan));

iNan = isnan(mX);
mX(iNan) = interp1(idx(~iNan),mX(~iNan),idx(iNan));

iNan = isnan(mY);
mY(iNan) = interp1(idx(~iNan),mY(~iNan),idx(iNan));

warning on

% Filter the data
sX(~isnan(sX)) = butter_filt(sX(~isnan(sX)),1,1/50,'low'); 
sY(~isnan(sY)) = butter_filt(sY(~isnan(sY)),1,1/50,'low'); 

mX(~isnan(mX)) = butter_filt(mX(~isnan(mX)),1,1/50,'low'); 
mY(~isnan(mY)) = butter_filt(mY(~isnan(mY)),1,1/50,'low'); 

% Store filtered data
body.statX  = sX;
body.statY  = sY;
body.mouthX = mX;
body.mouthY = mY;


function h = plotData(pl,cFrame,iPlate,displayAll,body,origin,S)

% Offset (in ix) for text
offset = 7;

% Define current plate
iPlate_c = iPlate;

hold on

% Include all plates, if displayAll
if displayAll
    iPlate = 1:length(pl);
end

% Plot body points
h(1,1) = plot(body.mouthX(cFrame),body.mouthY(cFrame),'mo');
h(2,1) = plot(body.statX(cFrame),body.statY(cFrame),'m+');
h(3,1) = plot([body.mouthX(cFrame) body.statX(cFrame)],...
            [body.mouthY(cFrame) body.statY(cFrame)],'m-');
h(4,1) = plot(origin(1),origin(2),'mo');

% Loop through plates to be displayed
for i = 1:length(iPlate)
    
    % Add fields 'angleX' and 'angleY', if not present in 'pl'
    if ~isfield(pl(iPlate(i)),'angleX') || ...
            isempty(pl(iPlate(i)).angleX)
        pl(iPlate(i)).angleX= nan(length(pl(iPlate(i)).tipX),1);
        pl(iPlate(i)).angleY= nan(length(pl(iPlate(i)).tipX),1);
    end

    % Define points for current plate
    bX = pl(iPlate(i)).baseX(cFrame);
    bY = pl(iPlate(i)).baseY(cFrame);
    tX = pl(iPlate(i)).tipX(cFrame);
    tY = pl(iPlate(i)).tipY(cFrame);
    aX = pl(iPlate(i)).angleX(cFrame);
    aY = pl(iPlate(i)).angleY(cFrame);
    
    % Plot line btwn base and angle points
    tmp = plot([aX bX],[aY bY],'w-');
    
    % Add tmp to handle list, clear tmp
    h = [h;tmp];
    clear tmp
    
    % Put emphasis on the current plate
    if iPlate_c == iPlate(i)
        tmp = plot([bX tX],[bY tY],'r-',bX,bY,'og',tX,tY,'r+');
    else
        tmp = plot([bX tX],[bY tY],'w-',bX,bY,'ow',tX,tY,'w+');
    end
    
    % Add tmp to handle list, clear tmp
    h = [h;tmp]; 
    clear tmp
    
    % Plot base point, if present
    if ~isnan(bX)
        
        tmp = text(bX-offset,bY+offset,...
                   num2str(pl(iPlate(i)).plate_num));
               
        % Add tmp to handle list, clear tmp
        h = [h;tmp];
        clear tmp
        
        % Put emphasis on the current base
        if iPlate_c == iPlate(i)
            set(h(end),'Color','g')
        else
            set(h(end),'Color','w')
        end
        
        % Put more emphasis, if keyframe
        if pl(iPlate(i)).baseMod(cFrame)
            set(h(end),'FontWeight','bold')
            set(h(end),'FontAngle','oblique')
            xlabel('Base key frame');
        else
            xlabel(' ');
        end
    end  
    
    clear bX bY tX tY aX aY num
end

function [cFrame,iPlate,pl] = continue_plate(pl)

% Make plate list
for i = 1:length(pl)
    plNums{i} = num2str(pl(i).plate_num);
end

% Prompt for selection of cureent plate number
[s,v] = listdlg('PromptString','Select plate number',...
    'SelectionMode','single',...
    'ListString',plNums);
if isempty(s)
    return
end

for i=1:length(pl)
    if pl(i).plate_num==str2num(plNums{s})
        iPlate = i;
        break
    end
end

% Define curent fraem from end of data
cFrame = find(~isnan(pl(iPlate).tipX),1,'last');
if isempty(cFrame)
    cFrame = 1;
end

clear s v

% Add angleX and angleY, if not present
if ~isfield(pl(iPlate),'angleX') || ...
        isempty(pl(iPlate).angleX)
    pl(iPlate).angleX= nan(seq.numFrames,1);
    pl(iPlate).angleY= nan(seq.numFrames,1);
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


