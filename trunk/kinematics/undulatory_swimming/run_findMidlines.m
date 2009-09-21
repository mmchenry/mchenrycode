function run_findMidlines
% Executes findMidlines on a group of sequences, then visualizes and runs
% stast on the results.  

%% Parameters

numAxis  = 15; % Number of point evaluation for axial kinematics
cutfreq  = 30;   % Low-pass filter cutoff freq (Hz)
numTailPts = 15; % Minimum number of tail points used to calculate wave speed

visSeries = 0;
visData   = 0;

%% Steps to execute

runMeasureBodyLength = 0;
runAcquireMidline    = 0;
runVisualizeMidlines = 0;
runFuseDirectories   = 1;
runPostProcessing1   = 0;
runExcludeOutliers   = 1;
runPostProcessing2   = 0;
runLinePlots         = 1;
runStats             = 1;

%% Get directories

%dataRoot    = '/Volumes/Docs/Projects/head_swimming/kinematic_data';
%vidRoot     = '/Users/mmchenry/Documents/MATLAB/shiner_video';

% ds          = data_dirs;
% dsVid       = video_dirs;

vidRoot     = '/Volumes/Docs/Projects/head_swimming/kinematic_data';
dataRoot    = '/Volumes/Docs/Projects/head_swimming/kinematic_data';

ds    = second_dirs;
dsVid = second_dirs;

%% Measure body lengths

if runMeasureBodyLength
    for i = 1:length(ds)      
        for j = 1:length(ds(i).dir)
            if ~isempty(ds(i).dir{j})
                disp(['loading i = ' num2str(i) ' j = ' num2str(j)]);
                
                load([dataRoot filesep ds(i).dir{j} filesep 'seq_params'])
                files = dir([vidRoot filesep dsVid(i).dir{j} filesep '*.tif']);
                
                im = imread([vidRoot filesep dsVid(i).dir{j} filesep ...
                             files(end).name]);
                disp('Pick off points along the midline');
                [x,y]   = choosePoints(im);
                
                if isempty(x)
                    bLength.skipMidline = 1;
                else
                    bLength.skipMidline = 0;
                    bLength.pix = max(cumsum((diff(x).^2+diff(y).^2).^0.5));
                    bLength.units = p.units;
                    dLength.in_units = bLength.pix .* p.calConst;
                end
                
                save([dataRoot filesep ds(i).dir{j} filesep ...
                        'bodyLength'],'bLength');
            end
        end
    end
end

%% Acquire midline
if runAcquireMidline
    
for i = 1:length(ds)
    
    disp(' ')
    disp(['Working on ' num2str(i) ' of ' num2str(length(ds)) ' fish'])
        
    for j = 1:length(ds(i).dir)
               
        if ~isempty(ds(i).dir{j}) && ...
            isempty(dir([dataRoot filesep ds(i).dir{j} filesep ...
                         'midline_data.mat']))
            
            % Load body length
            load([dataRoot filesep ds(i).dir{j} filesep 'bodyLength'])
            
            if ~bLength.skipMidline
                 % Status 
                tic
                disp(' ')
                disp(['    Working on ' num2str(j) ' of ' ...
                       num2str(length(ds(i).dir)) ' sequences . . .'])
                   
                % Load d
                load([dataRoot filesep ds(i).dir{j} filesep 'coord_data_posterior'])
                
                % Load p
                load([dataRoot filesep ds(i).dir{j} filesep 'seq_params'])
                    
                roi = [p.col_min p.row_min];  

                % Define video directory
                vid_dir = [vidRoot filesep dsVid(i).dir{j}];
                
                % Run findMidlines
                mid = findMidlines(vid_dir,bLength.pix,d,roi);
                
                Save
                save([dataRoot filesep ds(i).dir{j} filesep ...
                    'midline_data'],'mid');
                
                % Clear 
                clear mid vid_dir d
                
                % Status 
                tm = toc;
                disp(['             . . . done (' num2str(round(tm/60)) ' min).'])
                clear tm
            
            % Inform if skipMidlines is flagged
            elseif bLength.skipMidline
                disp(' ')
                disp(['    Skipping ' num2str(j) ' of ' ...
                       num2str(length(ds(i).dir)) ' (skipMidline==1)'])
            end
        
        % Inform if data already acquired
        elseif ~isempty(dir([dataRoot filesep ds(i).dir{j} filesep ...
                         'midline_data.mat']))
            disp(' ')
            disp(['    Skipping ' num2str(j) ' of ' ...
                       num2str(length(ds(i).dir)) ' (data collected)'])
        end
    end
    
    % Inform when done a fish
    disp(' ');
    disp(['Done ' num2str(i) ' of ' num2str(length(ds)) ' fish'])
end

end

%% Visualize midlines
if runVisualizeMidlines
    
% Parameters
numCols = 5;
numRows = 8;

% Create figure
hF = figure;
    
for i =  1:length(ds)
    
    k = 1;
    
    for j = 1:length(ds(i).dir)

        if ~isempty(ds(i).dir{j}) && ...
                ~isempty(dir([dataRoot filesep ds(i).dir{j} filesep ...
                'midline_data.mat']))
            
            % Change figure name
            set(hF,'Name',['Fish ' num2str(i) '  Seq ' num2str(j)])
            
            % Load d
            load([dataRoot filesep ds(i).dir{j} filesep 'coord_data_posterior'])

            % Load p
            load([dataRoot filesep ds(i).dir{j} filesep 'seq_params'])

            % Load midline data (mid)
            load([dataRoot filesep ds(i).dir{j} filesep 'midline_data'])

            % Get tiff file list, determine frame numbers
            vid_dir = dsVid(i).dir{j};
            files = dir([vidRoot filesep vid_dir filesep '*.tif']);
            numPanels = numCols*numRows;
            frNums = round(linspace(1,length(mid),numPanels));

            for k = 1:length(frNums)
                warning off
                % Load tiff
                im = imread([vidRoot filesep vid_dir filesep ...
                    files(frNums(k)).name]);
                
                % Plot image and midline
                subplot(numRows,numCols,k)
                imshow(im)
                hold on
                plot(mid(frNums(k)).mid(:,1),mid(frNums(k)).mid(:,2),'r-')
                title(['Frame ' num2str(frNums(k))])
                hold off
                warning on
                
            end
        end
    end
                
    % Prompt for return
    if k~=1
        disp(' ');disp('Press return to proceed')
        pause
    end
end

end

%% Fuse data directories
% Run to incorporate new sequences into the analysis
if runFuseDirectories
    clear ds  
    dRoot1 = data_dirs;
    dRoot2 = second_dirs;
    
    for i = 1:length(dRoot2)
        
        for j = 1:length(dRoot2(i).dir)
            if ~isempty(dRoot2(i).dir{j})
                dRoot1(i).dir{j} = dRoot2(i).dir{j};
            end
        end
        
    end
    
    ds = dRoot1;
    
    clear dRoot1
end

%% Post-processing 1: finding peaks
if runPostProcessing1
   
% Parameters
tolerance  = 3.e-1;
                
for i = 1:length(ds)
    
    disp(' ')
    disp(['Working on ' num2str(i) ' of ' num2str(length(ds)) ' fish'])
    
    k = 1;
    s = linspace(.001,.9,numAxis)';
    
    for j = 1:length(ds(i).dir)

        if ~isempty(ds(i).dir{j}) && ...
                ~isempty(dir([dataRoot filesep ds(i).dir{j} filesep ...
                'midline_data.mat'])) %&& ...
%                 isempty(dir([dataRoot filesep ds(i).dir{j} filesep ...
%                 'peak_data.mat']))
            
            % Directories
            vid_dir = [vidRoot filesep dsVid(i).dir{j}];
            imFiles = dir([vid_dir filesep '*.tif']);

            % Load p
            load([dataRoot filesep ds(i).dir{j} filesep 'seq_params'])
            framerate = p.framerate;
            
            clear p
            
            % Load midline data (mid)
            load([dataRoot filesep ds(i).dir{j} filesep 'midline_data'])         
            
            % Determine max size
            maxVals = 1;
            for k = 1:length(mid)
                maxVals = max([maxVals length(mid(k).s)]);
            end
            
            % Initiate matrices
            xVals = NaN(maxVals,length(mid));
            yVals = NaN(maxVals,length(mid));
            sVals = NaN(maxVals,length(mid));
            fVals = NaN(maxVals,length(mid));
            
            % Store values, stepping through frames
            for k = 1:length(mid)
                %xOrigin = mid(k).mid(1,1);
                xVals(1:length(mid(k).mid(:,1)),k) = mid(k).mid(:,1);
                yVals(1:length(mid(k).mid(:,2)),k) = mid(k).mid(:,2);
                sVals(1:length(mid(k).s),k)        = mid(k).s;
                fVals(:,k)                         = k.*ones(size(fVals,1),1);
            end
            clear k
            
            % Trim points, if  more than half the values are nans
            for k = 1:size(xVals,1)
                if sum(isnan(xVals(k,:))) > size(xVals,2)/2
                    xVals = xVals(1:k-1,:);
                    yVals = yVals(1:k-1,:);
                    sVals = sVals(1:k-1,:);
                    fVals = fVals(1:k-1,:);
                    break
                end
            end
            clear k
            
            % For each point, replace NaNs across time with mean of other values
            for k = 1:size(xVals,1)
                xVals(k,isnan(xVals(k,:))) = ...
                    mean(xVals(k,~isnan(xVals(k,:))));
                yVals(k,isnan(yVals(k,:))) = ...
                    mean(yVals(k,~isnan(yVals(k,:))));
                sVals(k,isnan(sVals(k,:))) = ...
                    mean(sVals(k,~isnan(sVals(k,:))));
            end
            clear k
            
            % Interpolate for consistant s
            sValsN = repmat(s,1,size(sVals,2));
            fValsN = repmat(fVals(1,:),length(s),1);
            xVals  = griddata(sVals,fVals,xVals,sValsN,fValsN);
            yVals  = griddata(sVals,fVals,yVals,sValsN,fValsN);
            
            sVals = sValsN;
            fVals = fValsN;
            
            clear sValsN fValsN
            
            % Check for NaNs
            if sum(isnan(xVals(:)))~=0
                % Repeat Trim points
                for k = 1:size(xVals,1)
                    if sum(isnan(xVals(k,:))) > size(xVals,2)/2
                        xVals = xVals(1:k-1,:);
                        yVals = yVals(1:k-1,:);
                        sVals = sVals(1:k-1,:);
                        fVals = fVals(1:k-1,:);
                        s     = s(1:k-1);
                        break
                    end
                end
                clear k
                
                % Repeat NaN replacement
                for k = 1:size(xVals,1)
                    xVals(k,isnan(xVals(k,:))) = ...
                        mean(xVals(k,~isnan(xVals(k,:))));
                    yVals(k,isnan(yVals(k,:))) = ...
                        mean(yVals(k,~isnan(yVals(k,:))));
                    sVals(k,isnan(sVals(k,:))) = ...
                        mean(sVals(k,~isnan(sVals(k,:))));
                end
                clear k
                
                % Confirm fix
                if sum(isnan(xVals(:)))~=0
                    error('NaN in xVals')
                end
            end
                
            % Step through each point and filter across time
            for k = 1:size(sVals,1)
                % Low-pass filter
                xVals_f(k,:) = butterworth(xVals(k,:),framerate,cutfreq,'low');
                yVals_f(k,:) = butterworth(yVals(k,:),framerate,cutfreq,'low');
                
                % High-pass filter
                xVals_f(k,:) = butter_high(xVals_f(k,:),framerate);
                yVals_f(k,:) = butter_high(yVals_f(k,:),framerate);
            end
            clear k
            
            % Surface smoothing splines for y coordinates
            spY  = spaps({s fVals(1,:)'},yVals_f,tolerance);
       
            % Body positions for evaluated splines 
            numSpPts = numAxis*20;
            sSp      = linspace(min(s),max(s),numSpPts);
            
            % Loop through time, find wave in the posterior region, store
            % in pk
            L = 1;
            M = 1;
            pk.s = [];
            for k = 1:size(xVals_f,2)
                
                % Evaluate spline for Y for current time
                Y = fnval(spY,{sSp k.*ones(size(sSp))});               
                Y = Y(:,1);
                
                % Fit cubic spline to current Y
                csY  = csapi(sSp,Y);
                
                % Find zeros for first dervative of Y. Interpolate for y
                % value at that peak
                zros  = fnzeros(fnder(csY,1));
                sWave = max(zros(:));
                yWave = interp1(sSp,Y,sWave);
                
                if ~isempty(zros)
                    % Locate coordinates of peak on video frame
                    xVid = interp1(sVals(:,k),xVals(:,k),sWave);
                    yVid = interp1(sVals(:,k),yVals(:,k),sWave);
                    
                    % Identify if new peak
                    if (k>1) && (~isempty(pk(end).s)) && ...
                            (sWave < pk(end).s(end))
                        L = L + 1;
                        M = 1;
                    end
                    
                    % Store position of peaks
                    pk(L).s(M,1)      = sWave;
                    pk(L).y(M,1)      = yWave;
                    pk(L).f(M,1)      = k;
                    pk(L).xVid(M,1)   = xVid;
                    pk(L).yVid(M,1)   = yVid;
                    pk(L).xMid_f(M,:) = xVals_f(:,k)';
                    pk(L).yMid_f(M,:) = yVals_f(:,k)';
                    pk(L).xMid(M,:)   = xVals(:,k)';
                    pk(L).yMid(M,:)   = yVals(:,k)';
                    pk(L).sMid(M,:)   = s';
                    
                    M = M + 1;
                    
                    if visSeries
                        subplot(2,1,1)
                        plot(sSp,Y,'k',sWave,yWave,'ro')
                        grid on
                        subplot(2,1,2)
                        im = imread([vid_dir filesep imFiles(k).name]);
                        warning off
                        imshow(im)
                        warning on
                        hold on
                        plot(xVid,yVid,'r+')
                        pause(.01)
                    end  
                end
                
                clear Y csY zros sWave yWave xVid yVid   
            end
            clear k L M
            
            % Plot results of finding peaks
            if visData
                for k = 1:length(pk)
                    if length(pk(k).s)>10
                        s = pk(k).s;
                        y = abs(pk(k).y);
                        plot(s,y,'o')
                        hold on
                        plot([min(s) max(s)],...
                            polyval(polyfit(s,y,2),[min(s) max(s)]),'k-')
                        xlabel('s'); ylabel('Y')
                    end
                end
                hold off
                title(['Fish ' num2str(i) ' seq ' num2str(j)])
                pause
            end
            clear k
            
            % Store results
            mid2.sVals = sVals;
            mid2.fVals = fVals;
            mid2.xVals = xVals;
            mid2.yVals = yVals;
            mid2.xVals_f = xVals_f;
            mid2.yVals_f = yVals_f;
            mid2.spY     = spY;
            
            % Save data
            save([dataRoot filesep ds(i).dir{j} filesep ...
                    'peak_data'],'pk');
                
            save([dataRoot filesep ds(i).dir{j} filesep ...
                    'midline_processed'],'mid2');
                
            % Clear values
            clear xVals yVals sVals fVals mid vid_dir imFiles mid2 pk
            clear xVals_f yVals_f sVals_f fVals_f spY numSpPts sSp 
        end
    end
    
    % Inform when done a fish
    disp(' ');
    disp(['Done ' num2str(i) ' of ' num2str(length(ds)) ' fish'])
end

clear s

end

%% Exclude outliers

if runExcludeOutliers
    ds(2).dir{6} = [];
    ds(7).dir{6} = [];
    ds(4).dir{4} = [];
    ds(4).dir{1} = [];
    ds(7).dir{1} = [];
end

%% Post-processing 2: Midline parameters
if runPostProcessing2

L_nosham = 1;MM = 1;cDir{1}=[];
L_sham = 1;

cats.pre_slow = [];  cats.pre_med = [];  cats.pre_fast = [];
cats.post_slow = []; cats.post_med = []; cats.post_fast = [];

d.ws = cats;
d.lambda = cats;
d.amp_coef1 = cats;
d.amp_coef2 = cats;
d.amp_coef3 = cats;
d.tba = cats;

dSh.ws = cats;
dSh.lambda = cats;
dSh.amp_coef1 = cats;
dSh.amp_coef2 = cats;
dSh.amp_coef3 = cats;
dSh.tba = cats;

clear cats

for i = 1:length(ds)
    
    for j = 1:length(ds(i).dir)

        if ~isempty(ds(i).dir{j}) && ...
                ~isempty(dir([dataRoot filesep ds(i).dir{j} filesep ...
                'peak_data.mat']))

            % Load p
            load([dataRoot filesep ds(i).dir{j} filesep 'seq_params'])
            
            framerate = p.framerate;
            calconst  = p.calconst;
            
            % Determine cSpd
            if ~isfield(p,'motorset')
                cSpd = 0;
                warning('No motorsetting given, assuming U_i_n_f = 0');
            elseif p.motorset==35
                cSpd = 4.5; %cm/s
            elseif p.motorset==40
                cSpd = 11; %cm/s
            elseif p.motorset==45
                cSpd = 22; %cm/s
            end
            
            clear p
            
            % Load bodyLength
            load([dataRoot filesep ds(i).dir{j} filesep 'bodyLength'])            
            bodyLength = bLength.pix;          
            clear bLength
            
            % Load pk
            load([dataRoot filesep ds(i).dir{j} filesep 'peak_data'])
            
            % Load m
           % load([dataRoot filesep ds(i).dir{j} filesep 'analyzed_data_posterior'])
            
            % Load snData (head kinematics)
            tmp = dir([dataRoot filesep ds(i).dir{j} filesep 'periodic_data_posterior.mat']);
            if ~isempty(tmp)
                load([dataRoot filesep ds(i).dir{j} filesep 'periodic_data_posterior'])
            else
                load([dataRoot filesep ds(i).dir{j} filesep 'periodic_data'])              
            end
            if i==6 && j==4
                tbf = p.tbf;
                clear p
            else
                tbf = snData.tbf;
                clear snData
            end
            
            
            % Load midline data (mid2)
            load([dataRoot filesep ds(i).dir{j} filesep 'midline_processed'])
            sVals = mid2.sVals(:,1);
            yVals = mid2.yVals;
            frames = mid2.fVals(1,:)';
            
            clear mid2
            
            % Step though each undulatory wave to find ws & lambda
            for k = 1:length(pk)
                % Grab data
                s = pk(k).s;
                t = pk(k).f .* (1/framerate);
                
                % Remove points that are not consecutive
                if sum(diff(pk(k).f)>1)>0
                    iJump = find(diff(pk(k).f)>1,1,'last');
                    s = s(iJump+1:end);
                    t = t(iJump+1:end);
                    %clear iJump
                end
                
                % Remove points that are not in the posterior
                iPost = s>0.7;
                s = s(iPost);
                t = t(iPost);
                
                
                % Skip, if fewer than 15 points
                if length(s)<15
                    iSkip(k)  = 1;
                    ws(k)     = 9999;
                    ws_r2(k)  = 9999;
                    lambda(k) = 9999;
                    
                % If more than 15 points, calc and store wavespeed
                else
                    % Calculate wave speed (ws)
                    coef = polyfit(t,s,1);
                
                    iSkip(k)  = 0;
                    ws(k)     = coef(1);
                    ws_r2(k)  = rSquared(s,polyval(coef,t));
                    
                    if 0 && j==6
                        disp([num2str(i) '   ' num2str(coef(1))])
                        plot(t,s,'o',t,polyval(coef,t),'k-')
                    ttt=2;
                    end
                    % Calculate wave length
                    lambda(k) = ws(k)/tbf;
                end
            end
            
            clear t s y k k2
            
            % Envelope 1: Define tail beat periods
            frInterval = round((1./tbf)/(1./framerate));
            ampPeriod  = 1:frInterval:length(frames);
            
            % Envelope 2: Take y-amplitude for each s over tailbeat periods
            for k = 1:length(ampPeriod)-1
                idx = (frames >= ampPeriod(k)) & (frames < ampPeriod(k+1));
                ys = yVals(:,idx);
                amps(:,k) = range(ys,2)/2;
            end
            
            % Envelope 3: Take envelope as mean of range values, ...
            % divided by body length
            amps  = mean(amps,2)./bodyLength;
            
            % Use nonlinear curve fit
            ampCoef  = nlinfit(sVals,amps,@myfun,[max(amps(:)) 0.3 .01]);
            amp_r2   = rSquared(amps,myfun(ampCoef,sVals));
            
            clear sVals amps ampPeriod frames frInterval ys yVals
 
            % Visualize model with coordinates over time
            if visSeries && j==6
                for k = 1:length(pk)
                    if ~iSkip(k)
                        t = pk(k).f .* (1/framerate);
                        
                        for n = 1:length(t)
                            phS = pi;
                            %x = pk(k).xMid_f(n,:)';
                            sAll = pk(k).sMid(n,:)';
                            yAll = pk(k).yMid_f(n,:)';
                            yPred = myfun(ampCoef,sAll) .* ...
                                sin((2*pi/lambda(k)).*(sAll-ws(k).*t(n))+phS);
                            
                            % Plot
                            plot(sAll,yAll./bodyLength,'o')
                            hold on
                            plot(sAll,yPred,'r-')
                            plot(sAll,myfun(ampCoef,sAll),'g--')
                            plot(sAll,-myfun(ampCoef,sAll),'g--')
                            hold off
                            title([num2str(i) ', ' num2str(j) '  ws = ' num2str(ws(k)) ...
                                '  lambda = ' num2str(lambda(k))])
                            pause(.01)
                        end
                    end
                end
            end
            
            
            % Store: General info
            if sum(~iSkip)>0
                if ds(i).pre(j)
                    exp_mode = 'pre ';
                else
                    exp_mode = 'post';
                end
                
                if ds(i).sham
                    L = L_sham;
                    dT = dSh;
                else
                    L = L_nosham;
                    dT = d;
                end
                
                sampleSize = length(~iSkip);
                fishnum = ds(i).fishnum;
                cDir{L} = ds(i).dir{j};
                
                % Store: wavespeed
                val = mean(ws(~iSkip));
                
                dT.ws = storeInCat(dT.ws,val,sampleSize,fishnum,j);
                
                dT.ws.all.exp(L,:)   = exp_mode;
                dT.ws.all.indiv(L,1) = fishnum;
                dT.ws.all.spd(L,1)   = cSpd;
                dT.ws.all.val(L,1)   = val;
                dT.ws.all.n(L,1)     = sampleSize;
                
                % Store: lambda
                val = mean(lambda(~iSkip));
                
                dT.lambda = storeInCat(dT.lambda,val,sampleSize,fishnum,j);
                
                dT.lambda.all.exp(L,:)   = exp_mode;
                dT.lambda.all.indiv(L,1) = fishnum;
                dT.lambda.all.spd(L,1)   = cSpd;
                dT.lambda.all.val(L,1)   = val;
                dT.lambda.all.n(L,1)     = sampleSize;
                
                % Store: amp_coef1
                val = ampCoef(1);
                
                dT.amp_coef1 = storeInCat(dT.amp_coef1,val,sampleSize,fishnum,j);
                
                dT.amp_coef1.all.exp(L,:)   = exp_mode;
                dT.amp_coef1.all.indiv(L,1) = fishnum;
                dT.amp_coef1.all.spd(L,1)   = cSpd;
                dT.amp_coef1.all.val(L,1)   = val;
                dT.amp_coef1.all.n(L,1)     = sampleSize;
                
                % Store: amp_coef2
                val = ampCoef(2);
                
                dT.amp_coef2 = storeInCat(dT.amp_coef2,val,sampleSize,fishnum,j);
                
                dT.amp_coef2.all.exp(L,:)   = exp_mode;
                dT.amp_coef2.all.indiv(L,1) = fishnum;
                dT.amp_coef2.all.spd(L,1)   = cSpd;
                dT.amp_coef2.all.val(L,1)   = val;
                dT.amp_coef2.all.n(L,1)     = sampleSize;
                
                % Store: amp_coef3
                val = ampCoef(3);
                
                dT.amp_coef3 = storeInCat(dT.amp_coef3,val,sampleSize,fishnum,j);
                
                dT.amp_coef3.all.exp(L,:)   = exp_mode;
                dT.amp_coef3.all.indiv(L,1) = fishnum;
                dT.amp_coef3.all.spd(L,1)   = cSpd;
                dT.amp_coef3.all.val(L,1)   = val;
                dT.amp_coef3.all.n(L,1)     = sampleSize;
                
                % Store: tba
                val = myfun(ampCoef,1);
                
                dT.tba = storeInCat(dT.tba,val,sampleSize,fishnum,j);
                
                dT.tba.all.exp(L,:)   = exp_mode;
                dT.tba.all.indiv(L,1) = fishnum;
                dT.tba.all.spd(L,1)   = cSpd;
                dT.tba.all.val(L,1)   = val;
                dT.tba.all.n(L,1)     = sampleSize;
                
                if ds(i).sham
                    L_sham = L_sham + 1;
                    dSh = dT;
                else
                    L_nosham = L_nosham + 1;
                    d = dT;
                end
                
                clear L dT
                
            end
            clear ws lambda tbf ampCoef amp_r2 bodyLength t exp_mode iSkip

        end
        
        if ~strcmp(cDir{end},ds(i).dir{j})
            nDir{MM} = ds(i).dir{j};
            MM = MM+1;
        end
    end
      
end

save([dataRoot filesep 'pooled_midline_parameters'],'d');
save([dataRoot filesep 'pooled_midline_parameters_Sham'],'dSh');

end

%% Line plots
if runLinePlots
    % Load d
    load([dataRoot filesep 'pooled_midline_parameters'])
    
    figure;
    
    subplot(2,3,1)
    lineplot(d.ws)
    axis square
    xlabel('speed')
    ylabel('Wave speed')
    
    subplot(2,3,2)
    lineplot(d.lambda)
    axis square
    xlabel('speed')
    ylabel('Wave length')
    
    subplot(2,3,3)
    lineplot(d.tba)
    axis square
    xlabel('speed')
    ylabel('tba')
    
    subplot(2,3,4)
    lineplot(d.amp_coef1)
    axis square
    xlabel('speed')
    ylabel('amp_coef1')
    
    subplot(2,3,5)
    lineplot(d.amp_coef2)
    axis square
    xlabel('speed')
    ylabel('amp coef2')
    
    subplot(2,3,6)
    lineplot(d.amp_coef3)
    axis square
    xlabel('speed')
    ylabel('amp coef 3')
end

%% Stats

if runStats
    
    load([dataRoot filesep 'pooled_midline_parameters'])
    
    % lambda
    s     = d.lambda.all;
    sName = 'lambda';
    
    varnames= {['exp_' sName] ['spd_' sName] ['indiv_' sName]};
    [p,table,stats] = anovan(s.val,{s.exp s.spd s.indiv},'model','interaction',...
        'varnames',varnames,'continuous',2);
    
    % tba
    s     = d.tba.all;
    sName = 'tba';
    
    varnames= {['exp_' sName] ['spd_' sName] ['indiv_' sName]};
    [p,table,stats] = anovan(s.val,{s.exp s.spd s.indiv},'model','interaction',...
        'varnames',varnames,'continuous',2);
    
    % coef1
    s     = d.amp_coef1.all;
    sName = 'amp_coef1';
    
    varnames= {['exp_' sName] ['spd_' sName] ['indiv_' sName]};
    [p,table,stats] = anovan(s.val,{s.exp s.spd s.indiv},'model','interaction',...
        'varnames',varnames,'continuous',2);
    
    % coef2
    s     = d.amp_coef2.all;
    sName = 'amp_coef2';
    
    varnames= {['exp_' sName] ['spd_' sName] ['indiv_' sName]};
    [p,table,stats] = anovan(s.val,{s.exp s.spd s.indiv},'model','interaction',...
        'varnames',varnames,'continuous',2);
    
    % coef3
    s     = d.amp_coef3.all;
    sName = 'amp_coef3';
    
    varnames= {['exp_' sName] ['spd_' sName] ['indiv_' sName]};
    [p,table,stats] = anovan(s.val,{s.exp s.spd s.indiv},'model','interaction',...
        'varnames',varnames,'continuous',2);
end



function d = storeInCat(d,val,sampleSize,fishnum,j)
if j==1
    d.pre_slow = [d.pre_slow; fishnum val sampleSize];
elseif j==2
    d.pre_med = [d.pre_med; fishnum val sampleSize];
elseif j==3
    d.pre_fast = [d.pre_fast; fishnum val sampleSize];
end
if j==4
    d.post_slow = [d.post_slow; fishnum val sampleSize];
elseif j==5
    d.post_med = [d.post_med; fishnum val sampleSize];
elseif j==6
    d.post_fast = [d.post_fast; fishnum val sampleSize];
end

function lineplot(d)

% Offset 
spd_off = .15;

% Colors used
clr_pre = [0.4784    0.0627    0.8941];
clr_post = [0    0.4980         0];

% Calculate stats for pre
[mu,sigma,muI,sigmaI] = normfit(d.pre_slow(:,2));
preD = [4.5-spd_off mu mu-muI(1) muI(2)-mu];

[mu,sigma,muI,sigmaI] = normfit(d.pre_med(:,2));
preD = [preD; 11-spd_off mu mu-muI(1) muI(2)-mu];

[mu,sigma,muI,sigmaI] = normfit(d.pre_fast(:,2));
preD = [preD; 22-spd_off mu mu-muI(1) muI(2)-mu];

% Calculate stats for post
[mu,sigma,muI,sigmaI] = normfit(d.post_slow(:,2));
postD = [4.5+spd_off mu mu-muI(1) muI(2)-mu];

[mu,sigma,muI,sigmaI] = normfit(d.post_med(:,2));
postD = [postD; 11+spd_off mu mu-muI(1) muI(2)-mu];

[mu,sigma,muI,sigmaI] = normfit(d.post_fast(:,2));
postD = [postD; 22+spd_off mu mu-muI(1) muI(2)-mu];

clear d

% Plot

h1 = plot(preD(:,1),preD(:,2),'-k');
hold on
h2 = plot(postD(:,1),postD(:,2),'-k');
set(h1,'Color',clr_pre);
set(h2,'Color',clr_post);
legend('pre','post')

h3 = errorbar(preD(:,1),preD(:,2),preD(:,3),preD(:,4),'k.');
set(h3,'Color',clr_pre);

h4 = errorbar(postD(:,1),postD(:,2),postD(:,3),postD(:,4),'k.');
set(h4,'Color',clr_post);

h5 = plot(preD(:,1),preD(:,2),'ok');
h6 = plot(postD(:,1),postD(:,2),'ok');
set(h5,'Color',clr_pre);
set(h6,'Color',clr_post);

set(gca,'XTick',[4.5 11 22])
set(gca,'XLim',[2.5 24])

function r2 = rSquared(X,Y)
% rsquared(X,Y).  This finds the r-squared value that describes the goodness of fit
% for a correlation between column vectors X and Y
if size(X,2)>1 | size(Y,2)>1
	error('You need to use two column vectors!');
end

% Example data from Chp 17 of Zar
%X = [3 4 5 6 8 9 10 11 12 14 15 16 17]';
%Y = [1.4 1.5 2.2 2.4 3.1 3.2 3.2 3.9 4.1 4.7 4.5 5.2 5.0]';

sig_x2 = sum(X.^2) - sum(X)^2/length(X);
sig_y2 = sum(Y.^2) - sum(Y)^2/length(Y);
sig_xy = sum(X.*Y) - sum(X).*sum(Y)/length(X);

ss_tot = sig_y2;
ss_regress = sig_xy^2/sig_x2;

r2 = ss_regress/ss_tot;

function y = myfun(beta,x)
% Normalize data before using this with nlinfit
y = beta(1).*(x-beta(2)).^2 + beta(3);

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

function data_filtered = butter_high(data,sample_rate)    
% All frequency values are in Hz.
Fs = sample_rate;   % Sampling Frequency

Fstop = .3;  % Stopband Frequency
Fpass = 2;   % Passband Frequency
Astop = 60;  % Stopband Attenuation (dB)
Apass = 1;   % Passband Ripple (dB)

% Calculate the order from the parameters using BUTTORD.
[N,Fc] = buttord(Fpass/(Fs/2), Fstop/(Fs/2), Apass, Astop);

% Calculate the zpk values using the BUTTER function.
[B A] = butter(N, Fc, 'high');

data_filtered   = filtfilt(B,A,data); 

function ds = data_dirs
% Structure of directories for data files

subDir = 'NeoFish1';
fnum = 1;
ds(fnum).fishnum   = fnum;
ds(fnum).sham      = 0;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'ControlFish3' filesep 'ControlFish3_3'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'ControlFish3' filesep 'ControlFish3_5'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'ControlFish3' filesep 'ControlFish3_9'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoFish1_1'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoFish1_4'];
ds(fnum).pre(6)    = 0;
ds(fnum).dir{6}    = [subDir filesep 'NeoFish1_7'];

subDir = 'NeoFish2';
fnum = 2;
ds(fnum).fishnum   = fnum;
ds(fnum).sham      = 0;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'ControlFish2' filesep 'ControlFish2_1'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'ControlFish2' filesep 'ControlFish2_7'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'ControlFish2' filesep 'ControlFish2_11'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoFish2_1'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoFish2_5'];
ds(fnum).pre(6)    = 0;
ds(fnum).dir{6}    = [subDir filesep 'NeoFish2_8'];

subDir = 'NeoFish3';
fnum = 3;
ds(fnum).fishnum   = fnum;
ds(fnum).sham      = 0;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'NeoFish3_1'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'NeoFish3_8'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'NeoFish3_10'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoFish3_12'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoFish3_16'];
ds(fnum).pre(6)    = 0;
ds(fnum).dir{6}    = [subDir filesep 'NeoFish3_18'];

subDir = 'NeoFish4';
fnum = 4;
ds(fnum).fishnum       = fnum;
ds(fnum).sham          = 0;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'NeoFish4_3'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'NeoFish4_5'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'NeoFish4_7'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoFish4_11'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoFish4_15'];
ds(fnum).pre(6)    = 0;
ds(fnum).dir{6}    = [subDir filesep 'NeoFish4_17'];

subDir = 'NeoFish5';
fnum = 5;
ds(fnum).fishnum   = fnum;
ds(fnum).sham      = 0;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'NeoFish5_3'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'NeoFish5_7'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'NeoFish5_10'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoFish5_14'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoFish5_16'];
ds(fnum).pre(6)    = 0;
ds(fnum).dir{6}    = [subDir filesep 'NeoFish5_19'];

subDir = 'NeoFish6';
fnum = 6;
ds(fnum).fishnum       = fnum;
ds(fnum).sham          = 0;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'NeoFish6_2'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'NeoFish6_5'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'NeoFish6_7'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoFish6_11'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoFish6_13'];
ds(fnum).pre(6)    = 0;
ds(fnum).dir{6}    = [];

subDir = 'NeoFish7';
fnum = 7;
ds(fnum).fishnum   = fnum;
ds(fnum).sham      = 0;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'NeoFish7_2'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'NeoFish7_7'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'NeoFish7_10'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoFish7_11'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoFish7_15'];
ds(fnum).pre(6)    = 0;
%ds(fnum).post_dir{3}   = [subDir filesep 'NeoFish7_18'];
ds(fnum).dir{6}   = [];

subDir = 'NeoControlFish1';
fnum = 8;
ds(fnum).fishnum   = fnum;
ds(fnum).sham      = 1;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'NeoControlFish1_2_S0001'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'NeoControlFish1_6_S0001'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'NeoControlFish1_8_S0001'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoControlFish1_10_S0001'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoControlFish1_14_S0001'];
ds(fnum).pre(6)    = 0;
ds(fnum).dir{6}    = [subDir filesep 'NeoControlFish1_17_S0001'];

subDir = 'NeoControlFish2';
fnum = 9;
ds(fnum).fishnum   = fnum;
ds(fnum).sham      = 1;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'NeoControlFish2_1'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'NeoControlFish2_4'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'NeoControlFish2_9'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoControlFish2_12'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoControlFish2_13'];
ds(fnum).pre(6)    = 0;
ds(fnum).dir{6}    = [subDir filesep 'NeoControlFish2_17'];

subDir = 'NeoControlFish3';
fnum = 10;
ds(fnum).fishnum  = fnum;
ds(fnum).sham     = 1;
ds(fnum).pre(1)   = 1;
ds(fnum).dir{1}   = [subDir filesep 'NeoControlFish3_1'];
ds(fnum).pre(2)   = 1;
ds(fnum).dir{2}   = [subDir filesep 'NeoControlFish3_4'];
ds(fnum).pre(3)   = 1;
ds(fnum).dir{3}   = [subDir filesep 'NeoControlFish3_9'];
ds(fnum).pre(4)   = 0;
ds(fnum).dir{4}   = [subDir filesep 'NeoControlFish3_11'];
ds(fnum).pre(5)   = 0;
ds(fnum).dir{5}   = [subDir filesep 'NeoControlFish3_14'];
ds(fnum).pre(6)   = 0;
ds(fnum).dir{6}   = [subDir filesep 'NeoControlFish3_16'];

function dsVid = video_dirs
% Structure of directories for video files

fnum = 1;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'ControlFish' filesep 'Controlfish3_03' '_S0001'];
dsVid(fnum).dir{2}    = [ 'ControlFish' filesep 'Controlfish3_05' '_S0001'];
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).dir{3}    = [ 'ControlFish' filesep 'Controlfish3_09' '_S0001'];
dsVid(fnum).dir{4}   = [ 'NeoFish' filesep 'NeoFish1_1' '_S0001'];
dsVid(fnum).dir{5}   = [ 'NeoFish' filesep 'NeoFish1_4' '_S0001'];
dsVid(fnum).dir{6}   = [ 'NeoFish' filesep 'NeoFish1_7' '_S0001'];

fnum = 2;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'ControlFish' filesep 'Controlfish2_1' '_S0001'];
dsVid(fnum).dir{2}    = [ 'ControlFish' filesep 'Controlfish2_7' '_S0001'];
dsVid(fnum).dir{3}    = [ 'ControlFish' filesep 'Controlfish2_11' '_S0001'];
dsVid(fnum).dir{4}   = [ 'NeoFish' filesep 'NeoFish2_1' '_S0001'];
dsVid(fnum).dir{5}   = [ 'NeoFish' filesep 'NeoFish2_5' '_S0001'];
dsVid(fnum).dir{6}   = [ 'NeoFish' filesep 'NeoFish2_8' '_S0001'];

fnum = 3;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'NeoFish' filesep 'NeoFish3_1' '_S0001'];
dsVid(fnum).dir{2}    = [ 'NeoFish' filesep 'NeoFish3_8' '_S0001'];
dsVid(fnum).dir{3}    = [ 'NeoFish' filesep 'NeoFish3_10' '_S0001'];
dsVid(fnum).dir{4}   = [ 'NeoFish' filesep 'NeoFish3_12' '_S0001'];
dsVid(fnum).dir{5}   = [ 'NeoFish' filesep 'NeoFish3_16' '_S0001'];
dsVid(fnum).dir{6}   = [ 'NeoFish' filesep 'NeoFish3_18' '_S0001'];

fnum = 4;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'NeoFish' filesep 'NeoFish4_3' '_S0001'];
dsVid(fnum).dir{2}    = [ 'NeoFish' filesep 'NeoFish4_5' '_S0001'];
dsVid(fnum).dir{3}    = [ 'NeoFish' filesep 'NeoFish4_7' '_S0001'];
dsVid(fnum).dir{4}   = [ 'NeoFish' filesep 'NeoFish4_11' '_S0001'];
dsVid(fnum).dir{5}   = [ 'NeoFish' filesep 'NeoFish4_15' '_S0001'];
dsVid(fnum).dir{6}   = [ 'NeoFish' filesep 'NeoFish4_17' '_S0001'];

fnum = 5;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'NeoFish' filesep 'NeoFish5_3' '_S0001'];
dsVid(fnum).dir{2}    = [ 'NeoFish' filesep 'NeoFish5_7' '_S0001'];
dsVid(fnum).dir{3}    = [ 'NeoFish' filesep 'NeoFish5_10' '_S0001'];
dsVid(fnum).dir{4}    = [ 'NeoFish' filesep 'NeoFish5_14' '_S0001'];
dsVid(fnum).dir{5}    = [ 'NeoFish' filesep 'NeoFish5_16' '_S0001'];
dsVid(fnum).dir{6}    = [ 'NeoFish' filesep 'NeoFish5_19' '_S0001'];

fnum = 6;
dsVid(fnum).fishnum   = fnum;
dsVid(fnum).sham      = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'NeoFish' filesep 'NeoFish6_2' '_S0001'];
dsVid(fnum).dir{2}    = [ 'NeoFish' filesep 'NeoFish6_5' '_S0001'];
dsVid(fnum).dir{3}    = [ 'NeoFish' filesep 'NeoFish6_7' '_S0001'];
dsVid(fnum).dir{4}    = [ 'NeoFish' filesep 'NeoFish6_11' '_S0001'];
dsVid(fnum).dir{5}    = [ 'NeoFish' filesep 'NeoFish6_13' '_S0001'];
dsVid(fnum).dir{6}    = [];

fnum = 7;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'NeoFish' filesep 'NeoFish7_2' '_S0001'];
dsVid(fnum).dir{2}    = [ 'NeoFish' filesep 'NeoFish7_7' '_S0001'];
dsVid(fnum).dir{3}    = [ 'NeoFish' filesep 'NeoFish7_10' '_S0001'];
dsVid(fnum).dir{4}   = [ 'NeoFish' filesep 'NeoFish7_11' '_S0001'];
dsVid(fnum).dir{5}   = [ 'NeoFish' filesep 'NeoFish7_15' '_S0001'];
%dsVid(fnum).dir{3}   = [subDir filesep'NeoFish' filesep  'NeoFish7_18'];
dsVid(fnum).dir{6}   = [];

fnum = 8;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 1;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'NeoControlFish' filesep 'NeoControlFish1_2_S0001'];
dsVid(fnum).dir{2}    = [ 'NeoControlFish' filesep 'NeoControlFish1_6_S0001'];
dsVid(fnum).dir{3}    = [ 'NeoControlFish' filesep 'NeoControlFish1_8_S0001'];
dsVid(fnum).dir{4}   = [ 'NeoControlFish' filesep 'NeoControlFish1_10_S0001'];
dsVid(fnum).dir{5}   = [ 'NeoControlFish' filesep 'NeoControlFish1_14_S0001'];
dsVid(fnum).dir{6}   = [ 'NeoControlFish' filesep 'NeoControlFish1_17_S0001'];

fnum = 9;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 1;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'NeoControlFish' filesep 'NeoControlFish2_1' '_S0001'];
dsVid(fnum).dir{2}    = [ 'NeoControlFish' filesep 'NeoControlFish2_4' '_S0001'];
dsVid(fnum).dir{3}    = [ 'NeoControlFish' filesep 'NeoControlFish2_9' '_S0001'];
dsVid(fnum).dir{4}   = [ 'NeoControlFish' filesep 'NeoControlFish2_12' '_S0001'];
dsVid(fnum).dir{5}   = [ 'NeoControlFish' filesep 'NeoControlFish2_13' '_S0001'];
dsVid(fnum).dir{6}   = [ 'NeoControlFish' filesep 'NeoControlFish2_17' '_S0001'];

fnum = 10;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 1;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'NeoControlFish' filesep 'NeoControlFish3_1' '_S0001'];
dsVid(fnum).dir{2}    = [ 'NeoControlFish' filesep 'NeoControlFish3_4' '_S0001'];
dsVid(fnum).dir{3}    = [ 'NeoControlFish' filesep 'NeoControlFish3_9' '_S0001'];
dsVid(fnum).dir{4}   = [ 'NeoControlFish' filesep 'NeoControlFish3_11' '_S0001'];
dsVid(fnum).dir{5}   = [ 'NeoControlFish' filesep 'NeoControlFish3_14' '_S0001'];
dsVid(fnum).dir{6}   = [ 'NeoControlFish' filesep 'NeoControlFish3_16' '_S0001'];

function dsVid = second_dirs
% vidRoot should be where all directories such as "ControlFish", "NeoFish"
% and "Tales" are present
%vidRoot     = '/Documents/Krijn/Shiners/';


fnum = 1;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'Tales' filesep 'Controlfish3_1'];
dsVid(fnum).dir{2}    = [];
dsVid(fnum).dir{3}    = [];
dsVid(fnum).dir{4}   = [];
dsVid(fnum).dir{5}   = [ ];
dsVid(fnum).dir{6}   = [ ];

 
fnum = 2;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ ];
dsVid(fnum).dir{2}    = [ ];
dsVid(fnum).dir{3}    = [ 'Tales' filesep 'Controlfish2_10'];
dsVid(fnum).dir{4}   = [ 'Tales' filesep 'NeoFish2_3'];
dsVid(fnum).dir{5}   = [ ];
dsVid(fnum).dir{6}   = [ 'Tales' filesep 'NeoFish2_8'];

 
fnum = 3;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'Tales' filesep 'NeoFish3_3'];
dsVid(fnum).dir{2}    = [];
dsVid(fnum).dir{3}    = [ 'Tales' filesep 'NeoFish3_10'];
dsVid(fnum).dir{4}    = [];
dsVid(fnum).dir{5}    = [ ];
dsVid(fnum).dir{6}    = [ ];

 
fnum = 4;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ ];
dsVid(fnum).dir{2}    = [ ];
dsVid(fnum).dir{3}    = [ 'Tales' filesep 'NeoFish4_7'];
dsVid(fnum).dir{4}   = [ 'Tales' filesep 'NeoFish4_11'];
dsVid(fnum).dir{5}   = [ ];
dsVid(fnum).dir{6}   = [ 'Tales' filesep 'NeoFish4_17'];

 
fnum = 5;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'Tales' filesep 'NeoFish5_1']; % may not work
%dsVid(fnum).dir{1}    = [ ];
dsVid(fnum).dir{2}    = [];
dsVid(fnum).dir{3}    = [ 'Tales' filesep 'NeoFish5_10'];
dsVid(fnum).dir{4}    = [ ];
dsVid(fnum).dir{5}    = [];
dsVid(fnum).dir{6}    = [];

 
fnum = 6;
dsVid(fnum).fishnum   = fnum;
dsVid(fnum).sham      = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'Tales' filesep 'NeoFish6_2'];
dsVid(fnum).dir{2}    = [ ];
dsVid(fnum).dir{3}    = [ 'Tales' filesep 'NeoFish6_7'];
dsVid(fnum).dir{4}    = [ 'Tales' filesep 'NeoFish6_11'];
dsVid(fnum).dir{5}    = [ ];
dsVid(fnum).dir{6}    = [];

 
fnum = 7;
dsVid(fnum).fishnum   = fnum;
dsVid(fnum).sham      = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'Tales' filesep 'NeoFish7_2']; %may not work
%dsVid(fnum).dir{1}    = [ ];
dsVid(fnum).dir{2}    = [ ];
dsVid(fnum).dir{3}    = [ ];
dsVid(fnum).dir{4}    = [ ];
dsVid(fnum).dir{5}    = [ ];
dsVid(fnum).dir{6}    = [ 'Tales' filesep 'NeoFish7_18'];

 







