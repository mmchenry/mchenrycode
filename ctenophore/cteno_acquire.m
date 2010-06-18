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

        clear answer
        
        
    case 'Continue an old comb plate'
        
        if isempty(a) || isempty(pl)
            error('No plate data exist')
        end
          
        % Make plate list
        for i = 1:length(pl)
            plNums{i} = num2str(pl(i).plate_num);
        end
  
        % Prompt for selection
        [s,v] = listdlg('PromptString','Select plate number',...
                        'SelectionMode','single',...
                        'ListString',plNums);
        if isempty(s)
            return
        end
        
        iPlate = str2num(plNums{s});
        
        cFrame = find(~isnan(pl(iPlate).tipX),1,'last');
        if isempty(cFrame)
            cFrame = 1;
        end
        
        clear s v
        
        
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
        
        iPlate = str2num(plNums{s});
        
        % Confirm delete
        but2 = questdlg(['Are you sure you want to delete plate ' ...
                        num2str(iPlate) '?'],'Warning','Yes - delete',...
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

baseModeColor = [.8 .4 .4];
frameSkip = 5;

% Give instructions
disp(' '); disp(' ');
disp('Left mouse   -  picks point.');disp(' ');
disp('Right mouse  - removes last point.');disp(' ');
disp('space bar - advance frame');disp(' ');
disp('left arrow - back one frame'); disp(' ');
disp('z - zoom mode (return to exit)');disp(' ');
disp('b - select base mode');disp(' ')
disp('t - select tip mode');disp(' ')
disp('j - jump to frame number');disp(' ')
disp('c - toggle color mode'); disp(' ')
disp(['+ - skip ' num2str(frameSkip) ' frames forward']); disp(' ')
disp(['- - skip ' num2str(frameSkip) ' frames backward']); disp(' ')
disp('s - change interval for skipping forward and backward');disp(' ')
disp('Press return when done collecting.')
disp('Press esc to exit');
disp(' '); disp(' ');

% Make figure window
hF = figure;
set(gcf,'DoubleBuffer','on');

% Set initial parameter values

numLimit = 10^10;
tipMode = 1;
colorMode = 0;

im = imread([imPath filesep seq.fileNames{1}]); 
hIm = imshow(im);

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
    
    if tipMode
        title(['Frame ' num2str(cFrame) ' of ' num2str(seq.numFrames)...
            '  TIP MODE'])
        set(hF,'Color','w')
    else
        title(['Frame ' num2str(cFrame) ' of ' num2str(seq.numFrames)...
            '  BASE MODE'])
        set(hF,'Color',baseModeColor);
    end
    
    
    % Loop for current frame
    while 1 
        
        % Plot existing data
        hold on
        h = plotData(pl,cFrame,iPlate);
        hold off
        
        [x,y,but] = ginput(1);
        
        % If return
        if isempty(but)
            break
            
        % Left click
        elseif but==1 
            
            if tipMode
                pl(iPlate).tipX(cFrame) = x;
                pl(iPlate).tipY(cFrame) = y;
            else
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
                
            end
            
            % Advance frame
            %cFrame = min([cFrame+1 seq.numFrames]);
            break
            
        % Right click    
        elseif but==3 
            
            if tipMode
                pl(iPlate).tipX(cFrame) = nan;
                pl(iPlate).tipY(cFrame) = nan;
            else
                
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
            %TODO: preserve zoom beyond frame advance
            
        % If 't'
        elseif but==116
            tipMode = 1;
            break
            
        % If 'b'
        elseif but==98
            tipMode = 0;
            break
            
        % If 'c'
        elseif but==99
            colorMode = abs(colorMode - 1);
            break
            
        % If 'j'
        elseif but==106
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
            
        % If spacebar or right arrow    
        elseif (but==32) || (but==29)           
            cFrame = min([cFrame+1 seq.numFrames]);
            break
            
        % If left arrow
        elseif but== 28
            cFrame = max([cFrame-1 1]);
            break
            
        % If -
        elseif but==45
            cFrame = max([cFrame-frameSkip 1]);
            break
            
        % If +   
        elseif but==43 || but==61         
            cFrame = min([cFrame+frameSkip seq.numFrames]);
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


function h = plotData(pl,cFrame,iPlate)

% Offset (in ix) for text
offset = 7;

% number of priot point to display
numPts = 5;

% % Plot earlier tips
% tX = pl(iPlate).tipX;
% tY = pl(iPlate).tipY;
% 
% tX = tX(1:cFrame);
% tY = tY(1:cFrame);
% 
% tX = tX(~isnan(tX));
% tY = tY(~isnan(tX));
% 
% idx = max([length(tX)-numPts 1]):1:(length(tX)-1);
% 
% h = plot(tX(idx),tY(idx),'r+');
%alpha(h,0.5);
hold on

% Plot current plate
bX = pl(iPlate).baseX(cFrame);
bY = pl(iPlate).baseY(cFrame);
tX = pl(iPlate).tipX(cFrame);
tY = pl(iPlate).tipY(cFrame);

num = num2str(pl(iPlate).plate_num);

h = plot([bX tX],[bY tY],'w-',bX,bY,'og',tX,tY,'r+');
%alpha(h1,0.5)
%plot();

if ~isnan(bX)
    hT = text(bX-offset,bY+offset,num);
    set(hT,'Color','g')
    if pl(iPlate).baseMod(cFrame)
        set(hT,'FontWeight','bold')
        set(hT,'FontAngle','oblique')
        xlabel('Base key frame');
    else
        xlabel(' ');
    end
end


