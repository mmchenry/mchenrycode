function track_single(movDir,show_steps)
% Tracks the movements of a single animal over a long duration.  Each step
% in the acquisition and analysis saves a .mat file in the movie directory
% so that step does not need to be repeated.  Delete these files if you want 
% to repeat a step.  The files names are meanImage.tif, seq_params.mat,
% coord_data.mat


%% Parameter values

% Path to first tif frame in sequence

% fPath = ...
% '/Volumes/workgroup/m_file_library/unfinished_karla/test_frames/Movie 3 00001.tif';

% fPath = ...
% ['/Volumes/workgroup/m_file_library/unfinished_karla/'...
% 'test_movie_for_matt3/test_movie_for_matt3 0001.tif'];

%fPath = '/Volumes/Docs/Projects/Karlas/Movie 3/Movie 3 00001.tif';
if nargin<1
    movDir = '/Volumes/My Book/Behavior Expts/20min_larva12';
end
%fPath = '/Volumes/workgroup/m_file_library/kinematics/single_tracker/Movie 5/Movie 5 00001.tif';

% fPath is the path to the first frame in a tiff sequence to be analyzed

% Mean images
maxFrames    = 100; % Maximum number of frames to be included in mean images
frInterval   = 10; % Interval over which to formulate a mean image
invert       = 1; % Whether to invert pixel values

% Analyzing light intensity
int_cutoff   = .3;   % Cut off for normalized light intensity 
min_gap      = 200;   % Min gap in frames btwn dark/light periods
diff_thresh  = 0.2;   % Normalized light change threshold

% Acquisition parameters
if nargin < 2
    show_steps = 0; % Display the acquisition
end
threshSpd    = 0.9; % Threshold body lengths/s for a tail beat
t_bin        = 20;  % Number of seconds to bin results
scale_factor = 5; % Radius of body roi in body lengths (3 is reasonable)
max_missed   = 10; % Max number of consecutively missed frames before error
saveInterval = 200; % Number of frames between saves
Darea_thresh = 2; % Max percent change in area allowed (0.5 is a good value)
displ_thresh = 2;  % Max number of body lengths accepted for a displacement

%% Add file separator to movDir

movDir = [movDir filesep];


% Gather info on movie
if fileThere('mov_info.mat',movDir)
    load([movDir filesep 'mov_info'])
    fPath = [movDir filesep mov.fileNames(1).name];
    mov.dirPath = movDir;
else
    a     = dir([movDir filesep '*.tif']);
    fPath = [movDir filesep a(1).name];
    mov   = findMov(fPath); 
end
save([movDir filesep 'mov_info.mat'],'mov')

endFrame     = mov.numFrames;
startFrame   = 1;



%% Gather light intensity data

if ~fileThere('intensity_data.mat',movDir)

    % Create wait bar
    disp(' '); disp(['Calculating intensity . . .']);
    h = waitbar(0,...
        ['Frame ' num2str(startFrame)],...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

    % Loop thru frames and calc mean intensity
    k = 1;
    for i = startFrame:endFrame
        [im.cdata,im.colormap] = ...
            imread([mov.dirPath filesep mov.fileNames(i).name '.' mov.ext]);
        %im          = grabFrame_nosub(mov,i,invert);
        frames(k)	= i;
        pix(k)      = mean(im.cdata(:));

        h = waitbar((i-startFrame)/(endFrame-startFrame),h,...
            ['Measuring light levels: Frame ' num2str(i) ' of ' ...
            num2str(endFrame-startFrame)]);
        if getappdata(h,'canceling')
            close force
            break
        end
        k = k+1;
    end
    
    % Normalize intensity data
    pix_norm   = (pix-min(pix))./range(pix);
    
    % Close waitbar
    close force
    
    % Save
    save([mov.dirPath filesep 'intensity_data'],'frames','pix','pix_norm');
    
    % Clear
    clear k h i 

else % If data exist, load

    disp(' ');disp('Using existing intensity data . . .');
    load([movDir filesep 'intensity_data.mat'])

end


%% Analyze light intensity

if ~fileThere('intensity_analyzed.mat',movDir)   
    
    % Define dark/light periods
    dk_frames   = frames(pix_norm<int_cutoff);
    lt_frames   = frames(pix_norm>(1-int_cutoff));
    dk_pix      = pix(pix_norm<int_cutoff);
    dk_pix_norm = pix_norm(pix_norm<int_cutoff);
    lt_pix      = pix(pix_norm>(1-int_cutoff));
    lt_pix_norm = pix_norm(pix_norm>(1-int_cutoff));
    
    % Determine if lights on at the beginning
    if min(lt_frames) < min(dk_frames)
        lt_on = 1;
    else
        lt_on = 0;
    end

    % Identify transition frames by looping through frames to find changes
    % in intensity
    iFrame  = 1;
    trans   = zeros(size(frames));

    while 1      
        % If lights go on
        if diff(pix_norm(iFrame:iFrame+1))>diff_thresh
            trans(iFrame) = 1;
            iFrame = iFrame + min_gap;        
        % If lights go off
        elseif diff(pix_norm(iFrame:iFrame+1))<(-diff_thresh)
            trans(iFrame) = -1;
            iFrame = iFrame + min_gap;           
        % Otherwise, advance
        else
            iFrame = iFrame + 1;
        end       
        % Break, if impossible
        if iFrame+1 > length(frames)
            break
        end
    end
    
    % Add final (trans)ition at the end
    if max(abs(trans(:)))==0 || min(abs(trans(:)))==1
        error('No light transitions detected');
    end
    trans(end) = -trans(find(abs(trans)==1,1,'last'));
    
    
    % Use dark/light frames and transitions to define cycles
    iCycle = 1;
    lastTrans = 1;
    while 1
        % If light is on
        if lt_on
            % Identify next transition 
            nextTrans  = frames(find(trans==-1,1));
            
            % Break, if near end
            if isempty(nextTrans)
                break
            end
            
            % Overwite current transition with a zero
            trans(find(trans==-1,1)) = 0;
            
            % Store data, update lt_on
            currFrames = ...
                lt_frames((lt_frames>=lastTrans)&(lt_frames<=nextTrans));
            currPix = lt_pix((lt_frames>=lastTrans)&(lt_frames<=nextTrans));
            
            c(iCycle).frames = currFrames;
            c(iCycle).pix    = currPix;
            c(iCycle).light_on = 1;
            lt_on = 0;
            iCycle = iCycle +1;
            lastTrans = nextTrans;
            
            clear currFrames currPix
            
        % If light is off
        else
            % Identify next transition 
            nextTrans  = frames(find(trans==1,1));
            
            % Break, if near end
            if isempty(nextTrans)
                break
            end
            
            % Overwite current transition with a zero
            trans(find(trans==1,1)) = 0;
            
            % Store data, update lt_on
            currFrames = ...
                dk_frames((dk_frames>=lastTrans)&(dk_frames<=nextTrans));
            currPix = ...
                dk_pix((dk_frames>=lastTrans)&(dk_frames<=nextTrans));
            
            c(iCycle).frames = currFrames;
            c(iCycle).pix    = currPix;
            c(iCycle).light_on = 0;
            lt_on = 1;
            iCycle = iCycle +1;
            lastTrans = nextTrans;
            
            clear currFrames currPix
            
        end
    end
    
    % Save data
    save([movDir filesep 'intensity_analyzed'],'c');
    clear im i h diff_thresh dk_frames iCycle iFrame int_cutoff lastTrans
    clear lt_frames lt_on min_gap nextTrans offColor onColor trans xMax
    clear xMin

else % If data exist, load

    disp(' ');disp('Using existing analyzed intensity data . . .');
    load([movDir filesep 'intensity_analyzed.mat'])
end


% Visualize how the code is identifying dark/ligh periods
if 0
    % Visualize light data
    offColor = .1.*[1 1 1];
    onColor  = .9.*[1 1 1];

    figure;
    subplot(2,1,1)
    for i = 1:length(c)
        xMin = min(c(i).frames);
        xMax = max(c(i).frames);

        if c(i).light_on
            h = fill([xMin xMin xMax xMax],[0 1 1 0],onColor);

        else
            h = fill([xMin xMin xMax xMax],[0 1 1 0],offColor);
        end
        set(h,'LineStyle','none')
        hold on
    end

    h = plot(frames,pix_norm,'r-');
    set(gca,'Color',[.5 .5 .5])
    set(h,'LineWidth',2)
    xlabel('Frames')
    ylabel('Normalized light intensity')

    clear xMin xMax h offColor onColor
end



%% Build list of light and dark frames

lt_k = 1;
dk_k = 1;


for i = 1:length(c)
    frameStep = round(length(c(i).frames)/maxFrames);
    
    if frameStep < frInterval
        error(['frame interval must be more than frame step.'...
            'Reduce maxFrames']);
    end
    
    for j = 1:frameStep:length(c(i).frames)
        % Check that interval doesn't extend too far
        if (c(i).frames(j)+frInterval-1)>c(i).frames(end)
            tmp = [c(i).frames(j):c(i).frames(end)]';
        else
            tmp = c(i).frames(j)+[0:frInterval-1]';
        end
        
        if c(i).light_on
            ltFrames(lt_k).f  = tmp;
            lt_k              = lt_k + 1;
        else
            dkFrames(dk_k).f  = tmp;
            dk_k              = dk_k + 1;
        end
    end
end

clear frameStep currFrames frInterval lt_k dk_k tmp


%% Prompt to see if larva moves in the first interval

if ~isfield(c(1),'movement')
    fr1 = imread([movDir filesep mov.fileNames(c(1).frames(1)).name '.tif']);
    fr2 = imread([movDir filesep mov.fileNames(c(1).frames(end)).name '.tif']);

    subplot(1,2,1)
    imshow(fr1); title(['Interval 1:  Frame ' num2str(c(1).frames(1))])
    subplot(1,2,2)
    imshow(fr2); title(['Interval 1:  Frame ' num2str(c(1).frames(end))])

    answ = questdlg('Does larva move in the initial interval?','Question',...
        'Yes','No','Cancel','Yes');

    if strcmp(answ,'Cancel')
        return

    elseif strcmp(answ,'No')
        c(1).movement = 0;

    else
        c(1).movement = 1;
    end

    save([movDir filesep 'intensity_analyzed'],'c');
end

%% Prompt to see if larva moves in the second interval

if ~isfield(c(2),'movement') || isempty(c(2).movement)
    fr1 = imread([movDir filesep mov.fileNames(c(2).frames(1)).name '.tif']);
    fr2 = imread([movDir filesep mov.fileNames(c(2).frames(end)).name '.tif']);

    subplot(1,2,1)
    imshow(fr1); title(['Interval 2:  Frame ' num2str(c(2).frames(1))])
    subplot(1,2,2)
    imshow(fr2); title(['Interval 2:  Frame ' num2str(c(2).frames(end))])

    answ = questdlg('Does larva move in the second interval?','Question',...
        'Yes','No','Cancel','Yes');

    if strcmp(answ,'Cancel')
        return
        
    elseif strcmp(answ,'Yes')
        c(2).movement = 1;
        
    elseif strcmp(answ,'No')
        error('Do not bother analyzing the larvae that never move');
    end

    save([movDir filesep 'intensity_analyzed'],'c');
end


%% Create dark mean image

if ~fileThere('meanImage_dk.tif',movDir)   
    
    h = waitbar(0,...
            ['Frame ' num2str(dkFrames(1).f(1))],...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    [img.cdata,tmp] = ...
         imread([movDir filesep mov.fileNames(1).name '.' mov.ext]);
    im          = zeros(size(img.cdata));
    
    for i = 1:length(dkFrames)
        
        for j = 1:length(dkFrames(i).f)
            frNum = dkFrames(i).f(j);
            [img.cdata,tmp] = ...
               imread([movDir filesep mov.fileNames(frNum).name '.' mov.ext]);
            im3D(:,:,j)  = img.cdata;
            %im3D(:,:,j)  = imadjust(img.cdata,stretchlim(img.cdata),[]);
        end
        
        im      = im + double(min(im3D,[],3));
        
        clear im3D img frNum
        
        %Update status
        h = waitbar(i/length(dkFrames),h,...
            ['Dark image: ' num2str(i) ' of ' ...
            num2str(length(dkFrames))]);
        if getappdata(h,'canceling')
            close force
            break
        end
    end
    
    im = uint8(round(im./length(dkFrames)));
    imwrite(im,[movDir 'meanImage_dk.tif'],'tif','Compression','none');
    
    close force
    clear im frames_curr k pDone h
else
    disp(' ');disp('Using existing dark mean image . . .');
end

clear dkFrames


%% Create light mean image

if ~fileThere('meanImage_lt.tif',movDir) 

    h = waitbar(0,...
            ['Frame ' num2str(ltFrames(1).f(1))],...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    [img.cdata,tmp] = ...
         imread([movDir filesep mov.fileNames(1).name '.' mov.ext]);
    im          = zeros(size(img.cdata));
    
    for i = 1:length(ltFrames)
        
        for j = 1:length(ltFrames(i).f)
            frNum = ltFrames(i).f(j);
            [img.cdata,tmp] = ...
               imread([movDir filesep mov.fileNames(frNum).name '.' mov.ext]);
            im3D(:,:,j)  = img.cdata;
            %im3D(:,:,j)  = imadjust(img.cdata,stretchlim(img.cdata),[]);
        end
        
        im      = im + double(min(im3D,[],3));
        
        clear im3D img frNum
        
        %Update status
        h = waitbar(i/length(ltFrames),h,...
            ['Light image: ' num2str(i) ' of ' ...
            num2str(length(ltFrames))]);
        if getappdata(h,'canceling')
            close force
            break
        end
    end
    
    im = uint8(round(im./length(ltFrames)));
    imwrite(im,[movDir 'meanImage_lt.tif'],'tif','Compression','none');
    
    close force
    clear im frames_curr k pDone h
else
    disp(' ');disp('Using existing light mean image . . .');
end

clear ltFrames


%% Select starting point, roi, threshold

if ~fileThere('seq_params.mat',movDir)

    yes_okay = 0;
    
    warning off all

    % Prompt for parameters
    prompt={'Frame rate (per sec)'};
    name='Parameters';
    numlines=1;
    defaultanswer={'13.36'};
    answer      = inputdlg(prompt,name,numlines,defaultanswer);
    if isempty(answer)
        return
    end
    p.framerate = str2num(answer{1});

    % Get coordinates for larva
    if c(1).light_on
        p.lt_on = 1;
    else
        p.lt_on = 0;
    end
    %img         = grabFrame(mov,c(1).frames(1),invert,p.lt_on);
    [img.cdata,img.colormap] = ...
        imread([movDir filesep mov.fileNames(1).name '.' mov.ext]);
    figure;
    [p.x,p.y] = choosePoints(img,0,'Select position of larva');
    close

    % Select dimensions of circular roi
    txt = 'Select vertical axis of roi';
    figure;
    [p.roi_v.x,p.roi_v.y]   = choosePoints(img,1,txt);

    txt = 'Select horizontal axis of roi';
    [p.roi_h.x,p.roi_h.y]   = choosePoints(img,1,txt);

    % Measure body length
    txt = 'Measure the body length';
    [xT,yT]   = choosePoints(img,1,txt);
    p.bLength = ((xT(2)-xT(1))^2 + (yT(2)-yT(1))^2)^0.5;
    clear xT yT txt
    close
    
    while ~yes_okay 
        
        % Grab frames for threshold finding
        if c(1).light_on
            img_lt  = grabFrame(mov,c(1).frames(1),invert,1);
            img_dk  = grabFrame(mov,c(2).frames(1),invert,0);
        else
            img_dk  = grabFrame(mov,c(1).frames(1),invert,0);
            img_lt  = grabFrame(mov,c(2).frames(1),invert,1);
        end
        
        % Matlab guesses a threshold value
        p.tVal_dk = graythresh(img_dk.cdata);
        p.tVal_lt = graythresh(img_lt.cdata);
        
        % Store path info in p
        p.path   = movDir;
        p.fname  = mov.fileNames(startFrame).name;
        
        % Run threshFinder to find threshold values
        % note: threshFinder saves p in seq_params.mat
        p.lt_on = 1;
        waitfor(threshFinder(img_lt.cdata,p))
        load([p.path filesep 'seq_params.mat'])
        
        p.lt_on = 0;
        waitfor(threshFinder(img_dk.cdata,p))
        load([p.path filesep 'seq_params.mat'])
         
        % Check that it's okay
        [x_roi,y_roi] = roiCoords(p);
        if c(1).light_on          
            c_on  = c(1);
            c_off = c(2);      
        else        
            c_on = c(2);
            c_off = c(1);
        end
        
        figure;
        
        % Light image (preview)
        im      = grabFrame(mov,c_on.frames(1),invert,1);
        imROI   = roipoly(im.cdata,x_roi,y_roi);
        imBW    = ~im2bw(im.cdata,im.colormap,p.tVal_lt);
        imBW    = imBW & imROI;
        subplot(1,2,1)
        imshow(imBW)
        title('Light interval:Sample frame');
        disp('You should be able to see only the larva');
        
        % Dark image (preview)
        im      = grabFrame(mov,c_off.frames(1),invert,0);
        imROI   = roipoly(im.cdata,x_roi,y_roi);
        imBW    = ~im2bw(im.cdata,im.colormap,p.tVal_dk);
        imBW    = imBW & imROI;
        subplot(1,2,2)
        imshow(imBW)      
        title('Dark interval:Sample frame');

        ButtonName = questdlg('Is these frames okay?', ...
            'Question', ...
            'Yes - proceed', 'No - redo', 'Cancel', 'Yes - proceed');
        switch ButtonName,
            case 'Yes - proceed',
                yes_okay = 1;
            case 'No - redo',
                yes_okay = 0;
            case 'Cancel',
                return

        end % switch
        
        close
        clear im imROI imBW x_roi y_roi c_on c_off
        warning on all

    end % while loop

else % if seq_param exists, load

    disp(' '); disp('Loading existing sequence parameters . . .'); 
    load([movDir filesep 'seq_params.mat'])

end

clear img



%% Acquire coordinates
%       Coodinate data is stored in d 

if fileThere('coord_data.mat',movDir)
    
    disp(' ');
    disp('Loading coordiates analyzed from coord_data.mat . . .')
    load([movDir filesep 'coord_data'])
  
else  % If run for the first time    
    
    d(1).x     = [];
    
end

% Skip acquisition if all frames analyzed
if (length(d)==length(c)) && (d(end).frame(end)==c(end).frames(end))
    
    disp(' ');
    disp(['All coordinates previously acquired. Skipping acquisition . . .'])

% Start acquisition, if frames remain unanalyzed
else  
    % Set startCycle (light cycle number at which to start acquisition)
    % & overwrite startFrame (frame num at which to initiate acquisition)
    currD = length(d);
    
    % If first cycle & no points yet
    if currD==1 && isempty(d(1).x)
        lastX    = p.x;
        lastY    = p.y;
        lastArea = 10^4; %Assumes larva will be largest blob
        
        if ~c(1).movement
            startCycle = 2;
        else
            startCycle = 1;
        end
        
        startFrame = c(startCycle).frames(1);
         
    %If not all points acquired in current cycle
    elseif length(d(currD).x) < length(c(currD).frames)
        startCycle = currD;
        startFrame = c(currD).frames(length(d(currD).x)+1);
        lastX      = d(currD).x(end);
        lastY      = d(currD).y(end);
        lastArea   = d(currD).area(end);
          
    % If all points are acquired, move onto the next cycle
    else
        startCycle = currD+1;
        startFrame = c(currD+1).frames(1);
        lastX      = d(currD).x(end);
        lastY      = d(currD).y(end);
        lastArea   = d(currD).area(end);
    end 
    clear currD
    
    % Set missedNum to zero
    missedNum  = 0; 
    successNum = 0;
    
    % Define elliptical roi
    [x_roi,y_roi] = roiCoords(p);
    
    % Make status bar
    if show_steps
        figure;
        set(gcf,'DoubleBuffer','on')
        pause(.3)
    else
        h = waitbar(0,...
            ['Frame ' num2str(startFrame)],...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    end
    
    
    % Loop through cycles__________________________________________________
    %      cNum is the light cycle number
    
    for cNum = startCycle:length(c)

        % Configure light_on parameter
        if c(cNum).light_on
            lt_on = 1;
            tVal    = p.tVal_lt;
        else
            lt_on = 0;
            tVal    = p.tVal_dk;
        end
        
        % Store cycle parameters
        d(cNum).lt_on = lt_on;

                
        % Loop through frames______________________________________________
        
        % Loop parameters
        %   If not all coords in current loop recorded, adjust index k
        if length(d(cNum).x) < length(c(cNum).frames)
            k = length(d(cNum).x)+1;
        else
            k          = 1;
        end
        
        for i = startFrame:c(cNum).frames(end)

            % Load frame, make correction, subtract outside dish roi
            im      = grabFrame(mov,i,invert,lt_on);
            %im.cdata= uint8(double(im.cdata)-correction);
            imROI   = roipoly(im.cdata,x_roi,y_roi);
            imBW    = ~im2bw(im.cdata,im.colormap,tVal);
            imBW    = imBW & imROI;
            clear imROI

            % Create circle roi around fish body
            [x_roiB,y_roiB] = roiBody(lastX,lastY,p.bLength,scale_factor);
            imROI   = roipoly(im.cdata,x_roiB,y_roiB);
            imBW2   = imBW & imROI;
            clear imROI

            % Dilate im & get properties
            se    = strel('disk',ceil(p.bLength/32),4);
            imBW3 = imdilate(imBW2,se);
            LL    = bwlabel(imBW3);
            props = regionprops(LL,'Centroid','Area');
            clear LL %imBW2           
            
            % If more than one point, choose closest in area ==============
            if length(props)>1 
                idx = 1;
%                 Darea_min = 10^10;
%                 dist_min  = 10^10;

                % Set min product value unrealistically high
                prod_min = 10^30;
                
                % Predicte next coordinates, if prior coordinates exist
                if isfield(d(cNum),'x') && (length(d(cNum).x)>2)
                    dX = d(cNum).x(end)-d(cNum).x(end-1);
                    dY = d(cNum).y(end)-d(cNum).y(end-1);
                    nextX = lastX + dX;
                    nextY = lastY + dY;
                else
                    nextX = lastX;
                    nextY = lastY;
                end
                %clear dX dY nextX nextY
                
                % Step through each blob
                for j=1:length(props)
                    %wght  = 2;
                    if k > 1
                        Darea = abs(props(j).Area - mean(d(cNum).area)); 
                    else
                        Darea = abs(props(j).Area - lastArea);   
                    end
                    
                    dist  = ((props(j).Centroid(1)-nextX)^2 + ...
                             (props(j).Centroid(2)-nextY)^2 )^0.5;
                         
                    prod  = Darea^3*dist;
                    
                    if prod < prod_min
                        idx = j;
                        %Darea_min = Darea;
                        prod_min = prod;
                        clear Darea
                    end
                    
                end
                
                if k > 1
                    Darea_min = abs(props(idx).Area - mean(d(cNum).area))...
                            ./mean(d(cNum).area);
                    displ = ((props(idx).Centroid(1)-lastX)^2 + ...
                             (props(idx).Centroid(2)-lastY)^2 )^0.5./p.bLength;
                end
                
                % Clear props (i.e. store no data), if area is above
                % threshold
                if (k > 1) && ...
                    ((Darea_min > Darea_thresh)) %|| (displ > displ_thresh))
                    props = [];
                    warning('Change in area too high')
                end
                
                clear Darea 
            
            % If only one point, store it =================================
            elseif length(props)==1 
                
                %Darea_min = abs((props(1).Area - lastArea))./lastArea;
                if k > 1
                    Darea_min = abs(props(1).Area - mean(d(cNum).area))...
                        ./mean(d(cNum).area);
                    displ = ((props(1).Centroid(1)-lastX)^2 + ...
                        (props(1).Centroid(2)-lastY)^2 )^0.5./p.bLength;
                end
                
                % Clear props (i.e. store no data), if area is above
                % threshold
                if (k > 1) && ...
                    ((Darea_min > Darea_thresh)) %|| (displ > displ_thresh))
                    props = [];
                    warning('Change in area too high')
                else
                    idx = 1;
                end
                
%                 % Store data
%                 d(cNum).x(k)     = props.Centroid(1);
%                 d(cNum).y(k)     = props.Centroid(2);
%                 d(cNum).area(k)  = props.Area;
%                 d(cNum).frame(k) = i;
%                 
%                 % Update for next loop
%                 missedNum  = 0;
%                 lastX      = d(cNum).x(k);
%                 lastY      = d(cNum).y(k);
%                 lastArea   = d(cNum).area(k);
%                 missedNum  = 0;
%                 successNum = successNum + 1;
%                 
%                 k          = k+1;
            end
            
            % If no points, warning =======================================
            if isempty(props) 
                
                warning(['Frame ' num2str(i) ' has no animal'])
                missedNum = missedNum+1;
                if missedNum > max_missed
                    beep; beep; beep;
                    error([num2str(max_missed) ...
                        ' consecutive frames have no animal']);
                end
                
            else
                
                % Store data
                d(cNum).x(k)     = props(idx).Centroid(1);
                d(cNum).y(k)     = props(idx).Centroid(2);
                d(cNum).area(k)  = props(idx).Area;
                d(cNum).frame(k) = i;
                %save([movDir filesep 'coord_data'],'d')
                
                % Update for next loop
                lastX      = d(cNum).x(k);
                lastY      = d(cNum).y(k);
                lastArea   = d(cNum).area(k);
                k          = k+1;
                missedNum  = 0;
                successNum = successNum + 1;
                
                clear idx
                
%                 if k>2
%                     disp([num2str(i) '  displ = ' num2str(displ) '   Darea = ' ...
%                         num2str(Darea_min)])
%                     clear Darea_min displ
%                 end
            end
            
                
                
            % Update status
            if show_steps %Update plot
                warning off all
                %imshow(imBW3)
                imshow(im.cdata)
                title(['Frame ' num2str(i)])
                hold on
                plot(lastX,lastY,'+r')
                plot(x_roiB,y_roiB,'b-')
                plot(x_roi,y_roi,'g--')
                hold off
                pause(0.01)
                %pause
                warning on all
            else% Update status bar
                h = waitbar((i)/(endFrame),h,...
                    ['Acquisition: Frame ' num2str(i) ' of ' ...
                    num2str(endFrame)]);
                if getappdata(h,'canceling')
                    
                    figure;
                    imshow(im.cdata)
                    title(['Frame ' num2str(i) ' (current)'])
                    hold on
                    plot(lastX,lastY,'+r')
                    plot(x_roiB,y_roiB,'b-')
                    
                    answ = questdlg('Resume acquisition?','Paused!',...
                        'Yes','No','Yes');
                    if strcmp(answ,'Yes')
                        close
                        setappdata(h,'canceling',0)
                        figure(h)
                    else
                        close
                        figure(h)
                        close force
                        break
                    end
                end
            end 
            
            % Save at regular interval
            if successNum >= saveInterval
                %disp('reenable save')
                save([movDir filesep 'coord_data'],'d')
                successNum = 0;
            end

        end % frame for-loop
        
        % Set startFrame for next cycle
        if cNum < length(c)
            startFrame = c(cNum+1).frames(1);
        end
        
    end % cycle for-loop
    

    % Close figure windows
    if show_steps
        close
    else
        close force
    end
    
    %Save when complete
    save([movDir filesep 'coord_data'],'d')
end


%% Analysis

if isempty(d(1).frame)
    startCycle = 2;
else
    startCycle = 1;
end

% Parameters
samplerate = p.framerate;
offColor = .3.*[1 1 1];
firstFrame = d(startCycle).frame(1);

figure;


% Loop through cycles ____________________________________

for i = startCycle:length(d)
    
    % Normalize position, define time
    t = (d(i).frame-firstFrame) ./ samplerate;
    x = d(i).x./p.bLength;
    y = d(i).x/p.bLength;

    % Calc rates
    t_spd   = t(2:end);
    spd     = ((diff(x).^2 + diff(y).^2).^0.5)./diff(t);
    t_accel = t_spd(1:end-1);
    accel   = diff(spd)./diff(t_spd);
    t_jerk  = t_accel(2:end);
    jerk    = diff(accel)./diff(t_accel);

    % Trimming first point(s) aligns spd, accel with jerks
    spd     = spd(3:end);
    accel   = accel(2:end);
    t       = t(3:end-1);

    % Find peaks
    idx    = find(spd > threshSpd & jerk > 0);
    spd_pk = spd(idx);
    t_pk   = t(idx);

    % Find peaks & speed within bins
    t_hist_tmp  = 0:t_bin:max(t);
    for j = 1:length(t_hist_tmp)-1
        idx1 = t>t_hist_tmp(j) & t<t_hist_tmp(j+1);
        if max(idx1)==0
            spd_bin(j)=0;
        else
            spd_bin(j) = mean(spd(idx1));
        end

        idx2 = t_pk>t_hist_tmp(j) & t_pk<t_hist_tmp(j+1);
        spikes(j) = sum(idx2);
        %if mean
    end

    % Shift t_hist_tmp for graphing
    t_hist  = t_hist_tmp(1:end-1) + diff(t_hist_tmp(1:2))/2;


    %Plot data
    subplot(3,1,1)
    plot(t./60,spd,'-')
    if ~c(i).light_on
        yLim = get(gca,'YLim');
        h = fill([min(t) min(t) max(t) max(t)]./60,...
             [yLim(1) yLim(2) yLim(2) yLim(1)],offColor);
        set(h,'LineStyle','none')
    end
    hold on
    plot(t./60,spd,'-',t_pk./60,spd_pk,'or')
    xlabel('time (min)')
    ylabel('speed (BL/s)');
    grid on
    xTicks = get(gca,'XTick');
    xLim   = get(gca,'XLim');

    subplot(3,1,2)
    if ~c(i).light_on
        bar(t_hist./60,spikes,1)
        yLim = get(gca,'YLim');
        h = fill([min(t) min(t) max(t) max(t)]./60,...
             [yLim(1) yLim(2) yLim(2) yLim(1)],offColor);
        set(h,'LineStyle','none')
    end
    hold on
    bar(t_hist./60,spikes,1)
    xlabel('time (min)');
    ylabel('number of spikes');
    grid on
    set(gca,'XTick',xTicks);
    set(gca,'XLim',xLim);

    subplot(3,1,3)
    if ~c(i).light_on
        bar(t_hist./60,spd_bin,1)
        yLim = get(gca,'YLim');
        h = fill([min(t) min(t) max(t) max(t)]./60,...
             [yLim(1) yLim(2) yLim(2) yLim(1)],offColor);
        set(h,'LineStyle','none')
    end
    hold on
    bar(t_hist./60,spd_bin,1)
    xlabel('time (min)');
    ylabel('mean speed (BL/s)');
    grid on
    set(gca,'XTick',xTicks);
    set(gca,'XLim',xLim);


end % Loop


% figure;
% plot(spikes./t_bin,spd_bin,'o')
% xlabel('Beat freq (Hz)')
% ylabel('Mean speed (BL/s)')











%% FUNCTIONS 

function data_filtered = butterworth(data,sample_rate,cutfreq,type)
%Returns data filtered by a butterworth filter
%  data - vector of data for filtering
%  sample_rate - rate of sampling for data (must have equivalent intervals)
%  cutfreq - cut off frequency (must be 2 values for bandstop)
%  type - 'low' for lowpass (default), 'high' for highpass, 'stop' for bandstop

if nargin < 4
    type = 'low';
end

if strcmp(type,'stop') && ~(length(cutfreq)==2)
    error('Stop pass must have a two-element cutfreq')
end 

ff              = cutfreq./(sample_rate./2);
[B A]           = butter(2,ff,type);
data_filtered   = filtfilt(B,A,data);     % Filtered data


function [x_roi,y_roi] = roiBody(x,y,bLength,scale_factor)
%Provides coordinates for an elliptical region of interest

numPts  = 200;
x_roi   = [];
y_roi   = [];
r       = bLength.*scale_factor;

theta   = linspace(0,pi/2,round(numPts/4))';
x_roi   = [x_roi; r .* cos(theta) + x];
y_roi   = [y_roi; r .* sin(theta) + y];

theta   = linspace(pi/2,pi,round(numPts/4))';
x_roi   = [x_roi; r .* cos(theta) + x];
y_roi   = [y_roi; r .* sin(theta) + y];

theta   = linspace(pi,1.5*pi,round(numPts/4))';
x_roi   = [x_roi; r .* cos(theta) + x];
y_roi   = [y_roi; r .* sin(theta) + y];

theta   = linspace(1.5*pi,2*pi,round(numPts/4))';
x_roi   = [x_roi; r .* cos(theta) + x];
y_roi   = [y_roi; r .* sin(theta) + y];


function [x_roi,y_roi] = roiCoords(p)
%Provides coordinates for an elliptical region of interest

numPts  = 400;
x_h     = p.roi_h.x(1:2);
y_v     = p.roi_v.y(1:2);
r_h     = abs(x_h(1)-x_h(2))/2;
r_v     = abs(y_v(1)-y_v(2))/2;
x_roi   = [];
y_roi   = [];

theta   = linspace(0,pi/2,round(numPts/4))';
x_roi   = [x_roi; r_h .* cos(theta) + mean(x_h)];
y_roi   = [y_roi; r_v .* sin(theta) + mean(y_v)];

theta   = linspace(pi/2,pi,round(numPts/4))';
x_roi   = [x_roi; r_h .* cos(theta) + mean(x_h)];
y_roi   = [y_roi; r_v .* sin(theta) + mean(y_v)];

theta   = linspace(pi,1.5*pi,round(numPts/4))';
x_roi   = [x_roi; r_h .* cos(theta) + mean(x_h)];
y_roi   = [y_roi; r_v .* sin(theta) + mean(y_v)];

theta   = linspace(1.5*pi,2*pi,round(numPts/4))';
x_roi   = [x_roi; r_h .* cos(theta) + mean(x_h)];
y_roi   = [y_roi; r_v .* sin(theta) + mean(y_v)];



function  y = fileThere(fName,fPath)
[tmp1,tmp2,ext,tmp3] = fileparts([fPath filesep fName]);
a	= dir([fPath filesep '*' ext]);
y	= 0;
for i = 1:length(a)
	if (length(a(i).name) > 3) && strcmp(a(i).name,fName)
		y = 1;
		break
	end
end


function img = grabFrame(mov,fNum,invert,lt_on)
%Return image data for a single video frame.  Accepts a TIFF stack or AVI.
%example: img = grabFrame(mov,fNum,invert)
%mov is the movie structure
%fNum is the frame number
%lt_on tells if lights are on

if mov.isTiff
    [img.cdata,img.colormap] = imread([mov.dirPath filesep mov.fileNames(fNum).name '.' mov.ext]);
else %assume AVI, if not TIFF
    img = (aviread([mov.dirPath filesep mov.fileName],fNum));
end

% if invert
%     img.colormap = img.colormap(end:-1:1,:);
% end

% Load subtraction image
if lt_on
    imSub  = imread([mov.dirPath filesep,'meanImage_lt.tif']);   
else
    imSub  = imread([mov.dirPath filesep,'meanImage_dk.tif']);
end

% Adjust grayscale values and convert to double
im     = img.cdata;
warning off
im     = imsubtract(imSub,im);
warning on

if invert
    im = imcomplement(im);
end

if min(im(:))==255
    warning(['Frame ' num2str(fNum) ' is completely white']);
else
    im = imadjust(im,[double(min(im(:)))/255;255/255],[10/255;255/255]);
end
img.cdata = im;


function mov = findMov(pName)
ext         = pName(max(find(pName=='.'))+1:end);
fName       = pName(max(find(pName==filesep))+1:end);
pName       = pName(1:max(find(pName==filesep)));
mov.dirPath = pName;
mov.ext     = ext;

%returns movie data 
if strcmp(ext,'avi')
    mov.isTiff      = 0;
    mov.fileName    = fName;
    mov.info        = aviinfo([mov.dirPath filesep mov.fileName]);
    mov.numFrames   = mov.info.NumFrames;
elseif strcmp(ext,'tif') || strcmp(ext,'tiff')
    mov.isTiff      = 1;
    mov.fileNames   = giveTiffStack(pName,fName);
    mov.numFrames   = length(mov.fileNames);
else
    error('Files should have either a .avi or .tif extension');
end


function files = giveTiffStack(mpath,fname)
% Returns a structure with info on the tiff stack.
% Assumes the last 5 digits are the file number

[pathstr,name,ext,versn]    = fileparts(fname);

% Determine the index iNum of the ending of the file name that 
% includes the frame number
%iZeros = find(name=='0');
%firstZero = min(find(name=='0'));
%iNum        = firstZero:length(name);
iNum =length(name)-4:length(name);

% define start of file name
%nameHead    = name(1:firstZero-1);
nameHead     = name(1:end-5);

% set up for loop
a           = dir(mpath);
startNum    = str2num(name(iNum:end));
tNum        = startNum;
j           = 1;
while 1==1
    nameEnd     = [num2str(zeros(1,length(iNum)-length(num2str(tNum)))) num2str(tNum)];
    tempName    = [nameHead nameEnd(find(~(nameEnd==' ')))];
    %tempName    = [nameHead nameEnd];
    isFile      = 0;
    %Step through file list to see if file exists:
    for i = (tNum-startNum)+1:length(a)
        [pathstr,name,ext,versn]    = fileparts(a(i).name);
        if ~(min(ext=='.')) & (length(name)>=length(tempName))
            if strcmp(name(1:length(tempName)),tempName)
                files(j).name   = tempName;
                j               = j + 1;
                tNum            = tNum + 1;
                isFile          = 1;
                break
            end
        end
    end
    if ~isFile, break; end
end

function [x,y] = choosePoints(img,link,txt)
%Used for finding coordinate points on a static image 'img'.
warning off all
imshow(img.cdata,img.colormap);
title(txt)
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
        imshow(img.cdata,img.colormap);
        title(txt)
        hold on
        if link
            plot(x,y,'ro-')
        else
            plot(x,y,'ro')
        end
    end
end


x = x'; y = y';
warning on all