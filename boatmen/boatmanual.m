function boatmanual(vPath)
% 
% vPath - path to tif files
% fName - filename of first frame
% 


%% Parameter values

% Header for image filenames
nameHead = 'seq';

% Extension for image files
nameSuffix = 'tif';

% Visualize acquisition
visSteps = 1;


num_digit = 6;

%% Get path of data file, load data

if nargin < 1
    vPath = uigetdir(pwd,'Select first frame');
end

% Load filenames for frames
a = dir([vPath filesep  '*' nameHead '*.' nameSuffix]);



%% Prompt for sequence info

% Look for mean image
a2 = dir([vPath filesep 'seq_info.mat']);

if isempty(a2)
    
    % Prompt for parameters
    prompt={'Frame rate (per sec)', ...
        'Start frame number',...
        'Last frame number'};
    name='Parameters';
    numlines=1;
    defaultanswer={'500','1',a(end).name(end-num_digit-length(nameSuffix):...
                            end-length(nameSuffix)-1),'',''};
    answer      = inputdlg(prompt,name,numlines,defaultanswer);
    if isempty(answer)
        return
    end
    
    p.framerate = str2num(answer{1});
    startFrame  = str2num(answer{2});
    endFrame    = str2num(answer{3}); 
    
    % Get indicies for video frames
    idx = 1;
    for i = 1:length(a)
        frNum = str2num(a(i).name(end-num_digit-length(nameSuffix):...
            end-length(nameSuffix)-1));     
        if (frNum >= startFrame) && (frNum <= endFrame) 
            p.frNums(idx) = frNum;
            p.filename{idx} = a(i).name;
            idx = idx + 1;
        end 
    end 

    warning on all

    save([vPath filesep 'seq_info.mat'],'p');
    
    clear startFrame endFrame prompt name numlines defaultanswer a
    
else % if seq_param exists, load

    disp(' '); disp('Loading existing starting point data . . .'); 
    load([vPath filesep 'seq_info.mat'])

end

clear img



%% Acquire landmarks from boatmen

% Color for highlighting 
%clr = [0.8 .6 0];
on_clr = [1 0 0];
off_clr = [1 1 1];

% Create figure window
f = figure;
set(f,'DoubleBuffer','on')
set(f,'CurrentCharacter','1')

% Check for data file
a3 = dir([vPath filesep 'boat_coords.mat']);

% Load py, if present
if ~isempty(a3)
    disp(' ');disp('Loading data . . .')
    load([vPath filesep 'boat_coords.mat'])
       
% Otherwise, create b
else
    
    im = imread([vPath filesep p.filename{1}]);
    warning off
    imshow(im)
    warning on
    
    b.xlim = xlim;
    b.ylim = ylim;
    
    % Fill data structure with nans
    nanfill = nan(length(p.frNums),1);
    
    b.xNose  = nanfill;
    b.yNose  = nanfill;
    b.xTail  = nanfill;
    b.yTail  = nanfill;
    b.xElbow = nanfill;
    b.yElbow = nanfill;
    b.xWrist = nanfill;
    b.yWrist = nanfill;
    b.xTip   = nanfill;
    b.yTip   = nanfill;
    
    clear nanfill
    
    % Overwrite nans with body points, if file present
    if ~isempty(dir([vPath filesep 'body_data.mat']))
        disp('Loading body_data.mat . . .')
        disp('')
        load([vPath filesep 'body_data.mat'])
        
        % Store previously acquired data
        b.xNose = body.headX;
        b.yNose = body.headY;
        b.xTail = body.tailX;
        b.yTail = body.tailY;
        
        clear body
    end
    
    % Overwrite nans with appendage points, if file present
    if ~isempty(dir([vPath filesep 'appendage_data.mat']))
        disp('Loading appendage_data.mat . . .')
        disp('')
        load([vPath filesep 'appendage_data.mat'])
        
        % Store previously acquired data
        b.xWrist  = pl(1).ptX;
        b.yWrist  = pl(1).ptY;
        
        clear body
    end
    
    % Current point capture mode
    b.cMode = 'elbow';
     
    % Current frame index
    if min(isnan(b.xNose))==1
        b.cIdx = 1;
    else
        b.cIdx = find(~isnan(b.xNose),1,'first');
    end

end


% Loop through frames, until finished
while 1
       
    % Display instructions
    disp(' ')
    disp('Commands:');
    disp('   Left arrow  - back up a frame')
    disp('   Right arrow - advance a frame')
    disp('   Right click - delete coordinate');
    disp('   Left click  - select new position for larva')
    disp('   "z"         - Zoom mode');
    disp('   "e"         - Elbow mode');
    disp('   "w"         - Wrist mode');
    disp('   "t"         - Tip mode');
    disp('   "n"         - Nose mode');
    disp('   "r"         - Rear mode');
    disp('   "j"         - Jump to another framez');
    disp(' ')
    
    while 1  
        
        % Update current frame number
        cFrame = p.frNums(b.cIdx);

        % Update image
        im = imread([vPath filesep p.filename{b.cIdx}]);
        
        % Display frame
        figure(f)   
        warning off all
        imshow(im)
        warning on all
        hold on
        title([b.cMode ' mode: Frame ' num2str(cFrame) ' of ' ...
            num2str(p.frNums(end)) ' (esc to quit)'])
        
        % Overlay body data
        hB = plot([b.xNose(b.cIdx) b.xTail(b.cIdx)],...
                  [b.yNose(b.cIdx) b.yTail(b.cIdx)],'-',...
                  b.xNose(b.cIdx),b.yNose(b.cIdx),'o');
             set(hB,'Color',off_clr);
             
        % Overlay elbow-wrist data
        h1 = plot([b.xElbow(b.cIdx) b.xWrist(b.cIdx)],...
                  [b.yElbow(b.cIdx) b.yWrist(b.cIdx)],'-');
             set(h1,'Color',off_clr);
             
        % Overlay tip-wrist data
        h2 = plot([b.xTip(b.cIdx) b.xWrist(b.cIdx)],...
                  [b.yTip(b.cIdx) b.yWrist(b.cIdx)],'-');
             set(h2,'Color',off_clr);
             
        % Highlight current mode
        if strcmp(b.cMode,'nose')
            set(hB(2),'Color',on_clr);
            
        elseif strcmp(b.cMode,'tail')
            hT = plot(b.xTail(b.cIdx),b.yTail(b.cIdx),'r+');
                 set(hT,'Color',on_clr)
                 
        elseif strcmp(b.cMode,'elbow')
            hT = plot(b.xElbow(b.cIdx),b.yElbow(b.cIdx),'r+');
                 set(hT,'Color',on_clr)
                 
        elseif strcmp(b.cMode,'wrist')
            hT = plot(b.xWrist(b.cIdx),b.yWrist(b.cIdx),'r+');
                 set(hT,'Color',on_clr)   
                 
        elseif strcmp(b.cMode,'tip')
            hT = plot(b.xTip(b.cIdx),b.yTip(b.cIdx),'r+');
                 set(hT,'Color',on_clr)   
                 
        end
             
        
        clear hT hB h1 h2
        
        hold off
        
        % Set limits (if zoomed)
        xlim(b.xlim);
        ylim(b.ylim);
        
        % Start coordinate acquire
        [cX,cY,cB] = ginput(1);
        
        % ESC
        if cB == 27
            close;
            return
            
%         % Return   
%         elseif isempty(cB)
%             close;
%             return
            
        % Left arrow
        elseif cB == 28
            % Back up frame
            if b.cIdx==1
                beep
            else
                b.cIdx = b.cIdx - 1;
            end
            
        % Right arrow
        elseif cB == 29
            % Advance frame
            if b.cIdx==length(p.frNums)
                beep;
            else
                b.cIdx = b.cIdx + 1;
            end
            
        % Right click
        elseif cB == 3 
            b = update_value(b,cX,cY,cB);
            
        % Left click
        elseif cB == 1        
            b = update_value(b,cX,cY,cB);
            
        % "z" -- zoom
        elseif cB == 122   
            title('Press return to exit zoom mode')
            zoom on
            pause
            b.xlim = xlim;
            b.ylim = ylim;        
            
        % "j" -- jump frame
        elseif cB == 106
            answer = questdlg('Jump to which frame?','','First','Last',...
                              'Frame number','First');
                          
            if isempty(answer)
                % Do nothing
                
            elseif strcmp(answer,'First')
                b.cIdx = 1;
                
            elseif strcmp(answer,'Last')
                b.cIdx = length(p.frNums);
                
            elseif strcmp(answer,'Frame number')
               answer = inputdlg({'Enter frame number'},'',1,...
                           {num2str(p.frNums(round(length(p.frNums)/2)))});        
               % Check input         
               if max(str2num(answer{1})==p.frNums)==0
                   warning('Requested frame number not an option');
                   
               else
                   b.cIdx = find(p.frNums==str2num(answer{1}),1,'first');
                   
               end
            end
            
        % "e" -- elbow mode
        elseif cB == 101      
            b.cMode = 'elbow';
            
            % "w" -- wrist mode
        elseif cB == 119    
            b.cMode = 'wrist';
            
            % "t" -- tip mode
        elseif cB == 116  
            b.cMode = 'tip';
            
            % "n" -- nose mode
        elseif cB == 110            
            b.cMode = 'nose';
            
            % "r" -- rear mode
        elseif cB == 114          
            b.cMode = 'rear';
            
        end

        clear im
        
        % Save data
        save([vPath filesep 'boat_coords'],'b')
        
    end
end

% Finished when capture reached
if (b.cFrame + 1 > length(p.frNums)) || ...
        ((p.captureFrame>0) && (p.frNums(b.cFrame)>=p.captureFrame))
    b.finished = 1;
end

% Advance index
b.cFrame = b.cFrame+1;

% Clear variables for next loop
clear img imBW imBW2 props imROI se x_roi y_roi maxB tmp






close
disp(' '); disp('     . . . Finished collecting data on boatman'); disp(' ')



function b = update_value(b,cX,cY,cB)

if strcmp(b.cMode,'nose')
    if cB==3
        b.xNose(b.cIdx)  = nan;
        b.yNose(b.cIdx)  = nan;      
    elseif cB==1
        b.xNose(b.cIdx)  = cX;
        b.yNose(b.cIdx)  = cY;
    end

elseif strcmp(b.cMode,'rear') 
    if cB==3
        b.xTail(b.cIdx)  = nan;
        b.yTail(b.cIdx)  = nan;
    elseif cB==1
        b.xTail(b.cIdx)  = cX;
        b.yTail(b.cIdx)  = cY;
    end


elseif strcmp(b.cMode,'elbow') 
    if cB==3
        b.xElbow(b.cIdx)  = nan;
        b.yElbow(b.cIdx)  = nan;
    elseif cB==1
        b.xElbow(b.cIdx)  = cX;
        b.yElbow(b.cIdx)  = cY;
    end
    
elseif strcmp(b.cMode,'wrist') 
    if cB==3
        b.xWrist(b.cIdx)  = nan;
        b.yWrist(b.cIdx)  = nan;
    elseif cB==1
        b.xWrist(b.cIdx)  = cX;
        b.yWrist(b.cIdx)  = cY;
    end    
    
elseif strcmp(b.cMode,'tip') 
    if cB==3
        b.xTip(b.cIdx)  = nan;
        b.yTip(b.cIdx)  = nan;
    elseif cB==1
        b.xTip(b.cIdx)  = cX;
        b.yTip(b.cIdx)  = cY;
    end     
    
end



