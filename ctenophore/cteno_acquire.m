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

% Look for data
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
    plNums{1} = ' ';
    
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



%% Acquisition mode

baseModeColor  = [.8 .4 .4];
angleModeColor = [.7 .7 .9];
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
disp( '  z - zoom mode (return to exit)'); 
disp( '  b - toggle base/tip mode');  
disp( '  c - toggle color mode'); 
disp( '  d - toggle displaying all comb plates'); 
disp( '  a - toggle angle mode');
disp( '  p - jump to another comb plate');
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

% Display image to get pixel coordinates
warning off
im  = imread([imPath filesep seq.fileNames{1}]); 
hIm = imshow(im);
warning on

% Get x and y limits
xlim_c = xlim;
ylim_c = ylim;


% Loop for interactive mode
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
    
    if progMode == 1
        title(['Frame ' num2str(cFrame) ' of ' num2str(seq.numFrames)...
            '  TIP MODE'])
        set(hF,'Color','w')
    elseif progMode == 0
        title(['Frame ' num2str(cFrame) ' of ' num2str(seq.numFrames)...
            '  BASE MODE'])
        set(hF,'Color',baseModeColor);
    elseif progMode == 2
        title(['Frame ' num2str(cFrame) ' of ' num2str(seq.numFrames)...
            '  ANGLE MODE'])
        set(hF,'Color',angleModeColor);
    end
    
    
    % Loop for current frame
    while 1 
        
        % Plot existing data
        hold on
        h = plotData(pl,cFrame,iPlate,displayAll);
        hold off
        
        [x,y,but] = ginput(1);
        
        % If return
        if isempty(but)
            break
            
        % Left click
        elseif but==1 
            
            if progMode == 1
                pl(iPlate).tipX(cFrame) = x;
                pl(iPlate).tipY(cFrame) = y;
                
            elseif progMode == 0
                % Set that the current base value has been modified
                pl(iPlate).baseMod(cFrame) = 1;
                
                % Locate slices after present that have been modified
                tmp = find(pl(iPlate).baseMod);
                idx = tmp((tmp>cFrame));
                
                % Define interval over which to define the base coords
                if isempty(idx)
                    iFill = min([cFrame seq.numFrames]):seq.numFrames;
                else
                    iFill = min([cFrame seq.numFrames]):idx(1);
                end
                
                % Update base coords
                pl(iPlate).baseX(iFill) = x.*ones(length(iFill),1);
                pl(iPlate).baseY(iFill) = y.*ones(length(iFill),1);
            
            elseif progMode == 2
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
                pl(iPlate).angleX(iFill) = x.*ones(length(iFill),1);
                pl(iPlate).angleY(iFill) = y.*ones(length(iFill),1);
                
            end
            
            % Advance frame
            %cFrame = min([cFrame+1 seq.numFrames]);
            break
            
        % Right click    
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
            elseif progMode == 2
                
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
                
            end
            
            % Back one frame
            %cFrame = max([cFrame-1 1]);
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
        elseif (but==98) || (but==116)
            if progMode == 1
                progMode = 0;
            elseif progMode == 0
                progMode = 1;
            end
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
            % Advance frame, if in base or tip mode
            if (progMode==1) || (progMode==0)
                cFrame = min([cFrame+1 seq.numFrames]);
                break
                
            % If angle mode, advace to next keyframe    
            elseif progMode==2
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
                
                break
            end
            
        % If left arrow
        elseif but== 28
            if (progMode==1) || (progMode==0)
            cFrame = max([cFrame-1 1]);
            break
            
            elseif progMode==2
                
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
                
                break
            end
            
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
            
        % If 'a'  
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
        
        % If 'p'
        elseif but==112
            [cFrame,iPlate,pl] = continue_plate(pl); 
            break
            
        end
        
        delete(h);
        
    end
    
    save([imPath filesep 'plate_data.mat'],'pl')
    
    if isempty(but) || (but==27)
        break
    end
    
end
close;


function h = plotData(pl,cFrame,iPlate,displayAll)

% Offset (in ix) for text
offset = 7;

% number of priot point to display
numPts = 5;

iPlate_c = iPlate;

hold on

% Include all plates, if displayAll
if displayAll
    iPlate = 1:length(pl);
end

h = [];

for i = 1:length(iPlate)
    
    % Add angleX and angleY, if not present
    if ~isfield(pl(iPlate(i)),'angleX') || ...
            isempty(pl(iPlate(i)).angleX)
        pl(iPlate(i)).angleX= nan(length(pl(iPlate(i)).tipX),1);
        pl(iPlate(i)).angleY= nan(length(pl(iPlate(i)).tipX),1);
    end

    % Plot current plate
    bX = pl(iPlate(i)).baseX(cFrame);
    bY = pl(iPlate(i)).baseY(cFrame);
    tX = pl(iPlate(i)).tipX(cFrame);
    tY = pl(iPlate(i)).tipY(cFrame);
    aX = pl(iPlate(i)).angleX(cFrame);
    aY = pl(iPlate(i)).angleY(cFrame);
    
    num = num2str(pl(iPlate(i)).plate_num);
    
    tmp = plot([aX bX],[aY bY],'w-');
    h = [h;tmp];
    clear tmp
    
    if iPlate_c == iPlate(i)
        tmp = plot([bX tX],[bY tY],'r-',bX,bY,'og',tX,tY,'r+');
    else
        tmp = plot([bX tX],[bY tY],'w-',bX,bY,'ow',tX,tY,'w+');
    end
    
    h = [h;tmp];
    clear tmp
    %alpha(h1,0.5)
    %plot();
    
    if ~isnan(bX)
        tmp = text(bX-offset,bY+offset,num);
        h = [h;tmp];
        clear tmp
        
        if iPlate_c == iPlate(i)
            set(h(end),'Color','g')
        else
            set(h(end),'Color','w')
        end
        
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