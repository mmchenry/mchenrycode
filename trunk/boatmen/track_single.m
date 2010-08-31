function track_single
% Tracks the movements of a single animal over a long duration.  Each step
% in the acquisition and analysis saves a .mat file in the movie directory
% so that step does not need to be repeated.  Delete these files if you want 
% to repeat a step.  The files names are meanImage.tif, seq_params.mat,
% coord_data.mat




%% Preliminaries

% Acquisition parameters

 fPath = ...
 '/Volumes/Lacie - K/Behavior Experiments/Movie 6/Movie 6 00001.tif';

% fPath is the path to the first frame in a tiff sequence to be analyzed

int_cutoff   = .3; % Cut off for normalized light intensity 
maxFrames   = 500; % Maximum number of frames to be included in mean images
invert       = 1; % Whether to invert pixel values
startFrame   = 1;
show_steps   = 1; % Display the acqusition
scale_factor = 3; % Radius of body roi in body lengths
max_missed   = 10; % Max number of consecutively missed frames before error

% Analyzing light intensity
min_gap     = 200;   % Min gap in frames btwn dark/light periods
diff_thresh = 0.2; % Normalized light change threshold

% Acquisition parameters
threshSpd  = 0.7; % Threshold body lengths/s for a tail beat
t_bin      = 4;  % Number of seconds to bin results

% Gather info on movie
mov      = findMov(fPath);
endFrame = mov.numFrames;



%% Gather light intensity data

if ~fileThere('intensity_data.mat',mov.dirPath)

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
    load([mov.dirPath filesep 'intensity_data.mat'])

end


%% Analyze light intensity

if ~fileThere('intensity_analyzed.mat',mov.dirPath)   
    
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
    
%     % Examine dark/light periods for drift, add a correction
%     c(1).pix_correction = 0;
%     c(2).pix_correction = 0;
%     if length(c)>2
%         for i = 3:2:length(c)
%             c(i).pix_correction = mean(c(i).pix)-mean(c(1).pix);
%         end
%     end
%     if length(c)>3
%         for i = 4:2:length(c)
%             c(i).pix_correction = mean(c(i).pix)-mean(c(2).pix);
%         end
%     end
    
    % Save data
    save([mov.dirPath filesep 'intensity_analyzed'],'c');
    clear im i h diff_thresh dk_frames iCycle iFrame int_cutoff lastTrans
    clear lt_frames lt_on min_gap nextTrans offColor onColor trans xMax
    clear xMin

else % If data exist, load

    disp(' ');disp('Using existing analyzed intensity data . . .');
    load([mov.dirPath filesep 'intensity_analyzed.mat'])
end

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


%% Create dark mean image

if ~fileThere('meanImage_dk.tif',mov.dirPath)    
    
    if ~c(1).light_on
        frames_curr = c(1).frames;
    else
        frames_curr = c(2).frames;
    end
    
%     [imCurr.cdata,imCurr.colormap] = ...
%             imread([mov.dirPath filesep mov.fileNames(frames_curr(1)).name '.' mov.ext]);
        
    imCurr      = grabFrame_nosub(mov,frames_curr(1),invert);
    im          = double(imCurr.cdata(:,:,1));
    
    disp(' '); disp(['Building dark mean image . . .']); 
    k = 10;
    
    for i = frames_curr(2):min([max(frames_curr) frames_curr(2)+maxFrames])
%         [imCurr.cdata,imCurr.colormap] = ...
%             imread([mov.dirPath filesep mov.fileNames(i).name '.' mov.ext]);
        
        imCurr  = grabFrame_nosub(mov,i,invert);
        im      = im + double(imCurr.cdata(:,:,1));
        pDone   = i./min([max(frames_curr) frames_curr(2)+maxFrames]) .* 100;
        if pDone > k
            disp(['  Done ' num2str(k) '%'])
            k = k + 10;
        end
    end
    
    disp(' '); disp('  All Done'); beep
    im = uint8(round(im./mov.numFrames));
    imwrite(im,[mov.dirPath 'meanImage_dk.tif'],'tif','Compression','none');
    
    clear im imCurr frames_curr k pDone
else
    disp(' ');disp('Using existing dark mean image . . .');
end


%% Create light mean image

if ~fileThere('meanImage_lt.tif',mov.dirPath)    
     
    if c(1).light_on
        frames_curr = c(1).frames;
    else
        frames_curr = c(2).frames;
    end
    
%     [imCurr.cdata,imCurr.colormap] = ...
%             imread([mov.dirPath filesep mov.fileNames(frames_curr(1)).name '.' mov.ext]);

    imCurr      = grabFrame_nosub(mov,frames_curr(1),invert);
    im          = double(imCurr.cdata(:,:,1));
    
    disp(' '); disp(['Building light mean image . . .']); 
    k = 10;
    
    for i = frames_curr(2):min([max(frames_curr) frames_curr(2)+maxFrames])
%         [imCurr.cdata,imCurr.colormap] = ...
%             imread([mov.dirPath filesep mov.fileNames(i).name '.' mov.ext]);
        imCurr  = grabFrame_nosub(mov,i,invert);
        im      = im + double(imCurr.cdata(:,:,1));
        pDone   = i./min([max(frames_curr) frames_curr(2)+maxFrames]) .* 100;
        if pDone > k
            disp(['  Done ' num2str(k) '%'])
            k = k + 10;
        end
    end
    
    disp(' '); disp('  All Done'); beep
    im = uint8(round(im./mov.numFrames));
    imwrite(im,[mov.dirPath 'meanImage_lt.tif'],'tif','Compression','none');
    
    clear im imCurr frames_curr k pDone
else
    disp(' ');disp('Using existing light mean image . . .');
end


clear maxFrames


%% Select starting point, roi, threshold

if ~fileThere('seq_params.mat',mov.dirPath)

    yes_okay = 0;

    while ~yes_okay
        warning off all
        
        % Prompt for parameters
        prompt={'Frame rate (per sec)'};
        name='Parameters';
        numlines=1;
        defaultanswer={'16.7'};
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
            imread([mov.dirPath filesep mov.fileNames(1).name '.' mov.ext]);
        figure;
        [p.x,p.y]   = choosePoints(img,0,'Select position of larva');
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
        p.path   = mov.dirPath;
        p.fname  = mov.fileNames(startFrame).name;
        
        % Run threshFinder to find threshold values
        % note: threshFinder saves p in seq_params.mat
        p.lt_on = 0;
        waitfor(threshFinder(img_dk.cdata,p))
        load([p.path filesep 'seq_params.mat'])
        
        p.lt_on = 1;
        waitfor(threshFinder(img_lt.cdata,p))
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
        title('Sample frame');
        disp('You should be able to see only the larva');
        
        % Dark image (preview)
        im      = grabFrame(mov,c_off.frames(1),invert,0);
        imROI   = roipoly(im.cdata,x_roi,y_roi);
        imBW    = ~im2bw(im.cdata,im.colormap,p.tVal_dk);
        imBW    = imBW & imROI;
        subplot(1,2,2)
        imshow(imBW)      

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
    load([mov.dirPath filesep 'seq_params.mat'])

end

clear img



%% Acquire coordinates

if fileThere('coord_data.mat',mov.dirPath)
    
    disp('Loading coordiates analyzed . . .')
    load([p.path filesep 'coord_data'])
  
else  % If run for the first time
    
    
    d(1).x     = [];
    
end

if (length(d)==length(c)) && (length(d(end).x)==endFrame)
    
    disp(['All coordinates previously acquired'])

else % Start acquisition, if frames remain
    
    % Set startCycle 
    currD = length(d);
    if length(d(currD).x) < length(c(currD).frames)
        startCycle = length(d);
    else
        startCycle = length(d)+1;  
    end 
    clear currD
    
    % Set lastX & lastY
    if startCycle==1
        lastX   = p.x;
        lastY   = p.y;
    else
        lastX   = d(startCycle-1).x(end);
        lastY   = d(startCycle-1).y(end);
    end
            
    % Set missedNum to zero
    missedNum = 0; 
    
    % Define elliptical roi
    [x_roi,y_roi] = roiCoords(p);
    
    % Make status bar
    if show_steps
        figure;
        set(gcf,'DoubleBuffer','on')
    else
        h = waitbar(0,...
            ['Frame ' num2str(startFrame)],...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    end
    
    
    % Loop through cycles__________________________________________________
    for cNum = startCycle:length(c)

        % Configure light_on parameter
        if c(startCycle).light_on
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
        cFrames    = c(cNum).frames;
        k          = 1;
        
        for i = cFrames

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
            se    = strel('disk',ceil(p.bLength/8),4);
            imBW3 = imdilate(imBW2,se);
            LL    = bwlabel(imBW3);
            props = regionprops(LL,'Centroid','Area');
            clear LL imBW2
            
            % Update status
            if show_steps %Update plot
                warning off all
                imshow(imBW)
                title(['Frame ' num2str(i)])
                hold on
                plot(lastX,lastY,'+r')
                plot(x_roiB,y_roiB,'b-')
                hold off
                pause(0.01)
                warning on all
            else% Update status bar
                h = waitbar((i-startFrame)/(endFrame-startFrame),h,...
                    ['Acquisition: Frame ' num2str(i) ' of ' ...
                    num2str(endFrame-startFrame)]);
                if getappdata(h,'canceling')
                    close force
                    break
                end
            end            
            
            if length(props)>1 % If more than one point, choose closest
                idx = 1;
                dist_min = 10^10;
                for j=1:length(props)
                    dist = ( (props(j).Centroid(1)-lastX)^2 +...
                        (props(j).Centroid(2)-lastY)^2 )^0.5;
                    if dist < dist_min
                        idx = j;
                        newX = props(j).Centroid(1);
                        newY = props(j).Centroid(2);
                        dist_min = dist;
                    end
                end
                
                % Store data
                d(cNum).x(k)     = props(idx).Centroid(1);
                d(cNum).y(k)     = props(idx).Centroid(2);
                d(cNum).frame(k) = i;
                save([p.path filesep 'coord_data'],'d')
                
                % Update for next loop
                lastX      = d(cNum).x(k);
                lastY      = d(cNum).y(k);
                k          = k+1;
                missedNum  = 0;

                clear idx dist newX newY dist_min

            elseif length(props)==1 % Store the one point in roi
                
                % Store data
                d(cNum).x(k)     = props.Centroid(1);
                d(cNum).y(k)     = props.Centroid(2);
                d(cNum).frame(k) = i;
                save([p.path filesep 'coord_data'],'d')
                
                % Update for next loop
                missedNum  = 0;
                lastX      = d(cNum).x(k);
                lastY      = d(cNum).y(k);
                k          = k+1;

            elseif isempty(props) % If no points
                
                warning(['Frame ' num2str(i) ' has no animal'])
                missedNum = missedNum+1;
                if missedNum > max_missed
                    beep; beep; beep;
                    error([num2str(max_missed) ...
                        ' consecutive frames have no animal']);
                end
                % Note: storing no data
            end

        end % frame for-loop
    end % cycle for-loop
    

    % Close figure windows
    if show_steps
        close
    else
        close force
    end


end


%% Analysis

% Parameters
samplerate = p.framerate;
offColor = .3.*[1 1 1];
firstFrame = d(1).frame(1);

figure;

% Loop through cycles ____________________________________
for i = 1:length(d)
    
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
        spd_bin(j) = mean(spd(idx1));

        idx2 = t_pk>t_hist_tmp(j) & t_pk<t_hist_tmp(j+1);
        spikes(j) = sum(idx2);
        %if mean
    end

    % Shift t_hist_tmp for graphing
    t_hist  = t_hist_tmp(1:end-1) + diff(t_hist_tmp(1:2))/2;


    %Plot data
    subplot(3,1,1)
    plot(t,spd,'-')
    if ~p(i).lt_on
        yLim = get(gca,'YLim');
        h = fill([min(t) min(t) max(t) max(t)],...
             [yLim(1) yLim(2) yLim(2) yLim(1)],offColor);
        set(h,'LineStyle','none')
    end
    hold on
    plot(t,spd,'-',t_pk,spd_pk,'or')
    xlabel('time (s)')
    ylabel('speed (BL/s)');
    grid on
    xTicks = get(gca,'XTick');
    xLim   = get(gca,'XLim');

    subplot(3,1,2)
    if ~p(i).lt_on
        bar(t_hist,spikes,1)
        yLim = get(gca,'YLim');
        h = fill([min(t) min(t) max(t) max(t)],...
             [yLim(1) yLim(2) yLim(2) yLim(1)],offColor);
        set(h,'LineStyle','none')
    end
    hold on
    bar(t_hist,spikes,1)
    xlabel('time (s)');
    ylabel('number of spikes');
    grid on
    set(gca,'XTick',xTicks);
    set(gca,'XLim',xLim);

    subplot(3,1,3)
    if ~p(i).lt_on
        bar(t_hist,spd_bin,1)
        yLim = get(gca,'YLim');
        h = fill([min(t) min(t) max(t) max(t)],...
             [yLim(1) yLim(2) yLim(2) yLim(1)],offColor);
        set(h,'LineStyle','none')
    end
    hold on
    bar(t_hist,spd_bin,1)
    xlabel('time (s)');
    ylabel('mean speed (BL/s)');
    grid on
    set(gca,'XTick',xTicks);
    set(gca,'XLim',xLim);


end % Loop


figure;
plot(spikes./t_bin,spd_bin,'o')
xlabel('Beat freq (Hz)')
ylabel('Mean speed (BL/s)')











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
a	= dir(fPath);
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
im     = (imadjust(img.cdata));
imSub  = (imadjust(imSub));

% Subtract image from current
% if size(img.cdata,3)>1
%     img.cdata(:,:,1)    = uint8(round(128 - imSub + im(:,:,1)));
%     img.cdata(:,:,2)    = uint8(round(128 - imSub  + im(:,:,2)));
%     img.cdata(:,:,3)    = uint8(round(128 - imSub  + im(:,:,3)));
% else
%     pix                 = 128 - imSub + im(:,:,1);
%     pix(find(pix>255))  = 255;
%     img.cdata(:,:,1)    = uint8(round(pix));
%     img.cdata
% end
warning off
im = imsubtract(imSub,im);
warning on

%im(find(im>255))  = 255;

if invert
    im = imcomplement(im);
end

img.cdata = im;

%img.colormap        = img.colormap(end:-1:1,:);



function img = grabFrame_nosub(mov,fNum,invert)
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

if invert
    img.colormap = img.colormap(end:-1:1,:);
end

% Adjust grayscale values
im     = imadjust(img.cdata);



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
elseif strcmp(ext,'tif') | strcmp(ext,'tiff')
    mov.isTiff      = 1;
    mov.fileNames   = giveTiffStack(pName,fName);
    mov.numFrames   = length(mov.fileNames);
else
    error('Files should have either a .avi or .tif extension');
end


function files = giveTiffStack(mpath,fname)
% Returns a structure with info on the tiff stack

[pathstr,name,ext,versn]    = fileparts(fname);

% Determine the index iNum of the ending of the file name that 
% includes the frame number
iZeros = find(name=='0');
% if max(diff(find(name=='0'))) > 1
%     firstZero = (iZeros(max(find(diff(iZeros)>1))+1));
% else
    firstZero = min(find(name=='0'));
%end
iNum        = firstZero:length(name);

% define start of file name
nameHead    = name(1:firstZero-1);

% set up for loop
a           = dir(mpath);
startNum    = str2num(name(iNum:end));
tNum        = startNum;
j           = 1;
while 1==1
    nameEnd     = [num2str(zeros(1,length(iNum)-length(num2str(tNum)))) num2str(tNum)];
    tempName    = [nameHead nameEnd(find(~(nameEnd==' ')))];
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