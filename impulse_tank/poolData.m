function poolData(mPath)

% mPath    - root path that holds the sequence directories (optional)
% 

%% Options for what to run

runCollect  = 1;
runPlotting = 1;
runStats    = 1;


%% Define directories
if ispc 
    
    mPath = '\\Flow.bio.uci.edu\workgroup\swimming_expts';
    
elseif nargin < 1
    
    % Define path for a mac with flow mounted
    if ~ispc && ~isempty(dir('/Volumes/workgroup/swimming_expts'))
        mPath = '/Volumes/workgroup/swimming_expts';
    else
        mPath = uigetdir(pwd,'Choose root path');
    end
    
end


%% Assemble batch list

if runCollect
    a = dir(mPath);
    
    % Loop thru directory, collect names of batch directories
    cBatch = [];
    k = 1;
    batch{1}=[];
    for i = 3:length(a)
        
        % Check if directory, starting with letter 'b'
        if a(i).isdir && length(a(i).name)> 6 && a(i).name(1)=='b'
            cBatch = [a(i).name(2:3) '/' a(i).name(4:5) '/' a(i).name(6:7)];
            
            % Loop through batch list to see if unique
            unique = 1;
            for j = 1:length(batch)
                if strcmp(cBatch,batch{j})
                    unique = 0;
                    break
                end
            end
            
            % Add to list, if unique
            if unique
                batch{k}.name = cBatch;
                batch{k}.filename = a(i).name;
                k = k+1;
            end
        end
    end
    
    if ~isfield(batch{1},'name')
        error('No batch directories within selected directory');
    end
end


%% Collect data from all sequences

if runCollect
    
    % Initiate k
    k = 1;
    
    % Initiate frequency matrix
    f = zeros(3,3,2);   
    
    % Initiate L, structure of latency values
    L.d_toward = [];
    L.d_side   = [];
    L.d_away   = [];
    L.a_beat   = [];
    L.a_glide  = [];
    L.a_still  = [];
    
    % Loop through batches
    for i = 1:length(batch)
        
        % Locate faststart_data.mat
        fileDir = dir([mPath filesep batch{i}.filename filesep 'faststart_data.mat']);
        
        % Give warning, skip batch, if data not in dir
        if isempty(fileDir)
            warning([mPath filesep batch{i}.filename ...
                filesep 'faststart_data.mat    not present']);
        
        % Otherwise, proceed
        else
            % Load fast start data (structure 'c')
            load([mPath filesep batch{i}.filename filesep 'faststart_data.mat'])
            
            % Loop through individuals in the batch
            for j = 1:length(c)
                
                % Check that data exist for larva
                if isempty(c(j).fast_start_frame_num) || ...
                   isempty(c(j).treshold_level) || ...
                   isempty(c(j).left_eye) || ...
                   isempty(c(j).right_eye) || ...
                   isempty(c(j).swimming) 
               
                    warning(['No data for larva ' num2str(j) ' in ' ...
                        batch{i}.filename]);
                    
                % If data exist, proceed    
                else
                    
                    % Extract fast start frame
                    faststart_frame = c(j).fast_start_frame_num;
                    
                    % Add flow_dir, if not there
                    if ~isfield(c(j),'flow_dir')
                        warning('flow_dir field missing, using R2L');
                        c(j).flow_dir = 'R2L';
                    end
                    
                    % Define vector for flow direction
                    if strcmp(c(j).flow_dir,'L2R')
                        flowVector = [1 0];
                    elseif strcmp(c(j).flow_dir,'R2L')
                        flowVector = [-1 0];
                    elseif strcmp(c(j).flow_dir,'D2U')
                        flowVector = [0 1];
                    elseif strcmp(c(j).flow_dir,'U2D')
                        flowVector = [0 -1];
                    end
                    
                    % Load coords for the eyes, calc vector
                    lEye = [c(j).left_eye.x c(j).left_eye.y];
                    rEye = [c(j).right_eye.x c(j).right_eye.y];
                    eyeVector = lEye - rEye;
                    clear lEye rEye
                    
                    % Calculate angle
                    theta = (180/pi) * acos(dot(eyeVector,flowVector)/...
                        norm(eyeVector)*norm(flowVector));
       %TODO: Change these categories             
                    if (theta > 0) && (theta <= 60)
                        direction = 'toward';
                    elseif (theta > 60) && (theta <= 120)
                        direction = 'side  ';
                    elseif (theta > 120) && (theta <= 180)
                        direction = 'away  ';
                    end
                    
                    % Calculate latency
                    [amp,daq,vid] = synchronizeData(...
                                        [mPath filesep batch{i}.filename],0);
                    preTrig_frms  = vid.dur - vid.postTrg;
                    
                    % Define start frame
                    start_frm = c(j).fast_start_frame_num;
                    
                    % If no fast start
                    if start_frm==0
                        
                        fast_start = 2;
                    
                    % If fast start, calculate latency
                    else
  
                        % Latency is the video time for the frame at which
                        % motion was visible - video time of pretrigger period
                        % - the delay in the satrt of the motor motion
                        latency = vid.t(start_frm+preTrig_frms) - amp.t_mStart;                        
                        
                        if latency < 0
                            warning(['start_frm = ' num2str(start_frm)...
                                ' in fish ' num2str(j) ' ' ...
                                batch{i}.filename ...
                                ' is before stimulus '...
                                '-- not recorded']);
                            fast_start = 0;
                        else
                            
                            % Store data for latency ANOVA
                            s.activity(k,:)   = c(j).swimming;
                            s.switchTime(k,1) = c(j).switchTime;
                            s.direction(k,:)  = direction;
                            s.latency(k,1)    = latency;
                            
                            k = k + 1;
                            
                            % Store data for latency plotting
                            if strcmp(direction,'toward')
                                L.d_toward = [L.d_toward; latency];
                            elseif strcmp(direction, 'side  ')
 
                                L.d_side = [L.d_side; latency];
                            elseif strcmp(direction, 'away  ')
                                L.d_away = [L.d_away; latency];
                            end
                            
                            if strcmp(c(j).swimming,'beats')
                                L.a_beat = [L.a_beat; latency];
                            elseif strcmp(c(j).swimming,'glide')
                                L.a_glide = [L.a_glide; latency];
                            elseif strcmp(c(j).swimming,'still')
                                L.a_still = [L.a_still; latency];
                            end    
                            
                            % Set fast_start
                            fast_start = 1;
                        end
                    end
                    
                    
                    % Record response (or non-response)
                    if (fast_start==1) || (fast_start==2)             
                        
                        if strcmp(c(j).swimming,'beats')
                            if strcmp(direction,'toward')
                                f(1,1,fast_start) = f(1,1,fast_start) + 1;
                            elseif strcmp(direction, 'side  ')
                                f(1,2,fast_start) = f(1,2,fast_start) + 1;
                            elseif strcmp(direction, 'away  ')
                                f(1,3,fast_start) = f(1,3,fast_start) + 1;
                            end
                        elseif strcmp(c(j).swimming,'glide')
                            if strcmp(direction,'toward')
                                f(2,1,fast_start) = f(2,1,fast_start) + 1;
                            elseif strcmp(direction, 'side  ')
                                f(2,2,fast_start) = f(2,2,fast_start) + 1;
                            elseif strcmp(direction, 'away  ')
                                f(2,3,fast_start) = f(2,3,fast_start) + 1;
                            end
                        elseif strcmp(c(j).swimming,'still')
                            if strcmp(direction,'toward')
                                f(3,1,fast_start) = f(3,1,fast_start) + 1;
                            elseif strcmp(direction, 'side  ')
                                f(3,2,fast_start) = f(3,2,fast_start) + 1;
                            elseif strcmp(direction, 'away  ')
                                f(3,3,fast_start) = f(3,3,fast_start) + 1;
                            end
                        end
                        
                    end                           
                    
                    clear amp daq vid fast_start
                    clear direction theta latency
                    
                end
            end            
        end
    end
    
    save([mPath filesep 'ANOVA_data.mat'],'s')
    save([mPath filesep 'probability_data.mat'],'f')
    save([mPath filesep 'latency_data.mat'],'L')
    
end


%% Visualize results

if runPlotting
    
    % Load L, latency data
    load([mPath filesep 'latency_data.mat'])
    
    % Load f
    load([mPath filesep 'probability_data.mat'])
    load([mPath filesep 'ANOVA_data.mat'])
    
    % Parameters for plotting
    %yTicksL = (floor(min(1000.*s.latency)):10:ceil(max(1000.*s.latency)))./1000;
    yTicksP = 0:.5:1;
    txtOffLatency = .004;
    txtOffProb = 0.05;
    
    % Calculate probability stats
    [pToward,iToward] = binofit(sum(f(:,1,1)),sum(sum(f(:,1,:))));
    iToward = iToward(2) - pToward;
    
    [pAway,iAway] = binofit(sum(f(:,2,1)),sum(sum(f(:,2,:))));
    iAway = iAway(2) - pAway;
    
    [pSide,iSide] = binofit(sum(f(:,3,1)),sum(sum(f(:,3,:))));
    iSide = iSide(2) - pSide;
    
    [pBeat,iBeat] = binofit(sum(f(1,:,1)),sum(sum(f(1,:,:))));
    iBeat = iBeat(2) - pBeat;
    
    [pGlide,iGlide] = binofit(sum(f(2,:,1)),sum(sum(f(2,:,:))));
    iGlide = iGlide(2) - pGlide;
    
    [pStill,iStill] = binofit(sum(f(3,:,1)),sum(sum(f(3,:,:))));
    iStill = iStill(2) - pStill;
    
    figure;
    
    % Plot probability vs. activity level
    subplot(2,2,1)
        % Plot
        errorbar([1 2 3],[pBeat; pGlide; pStill],[iBeat; iGlide; iStill],'+k') 
        hold on
        bar([1 2 3],[pBeat; pGlide; pStill],'w')
    
        % Labels & axes
        xlabel('beat                         glide                         still')
        set(gca,'XTickLabel',' ')
        set(gca,'YTick',yTicksP)
        set(gca,'YLim',[yTicksP(1) yTicksP(end)]);
        ylabel('Probability')
        hold off
        
        % Numbers
        text(1,pBeat-txtOffProb,[num2str(sum(f(1,:,1))) '/' ...
                                    num2str(sum(sum(f(1,:,:))))]);
        text(2,pGlide-txtOffProb,[num2str(sum(f(2,:,1))) '/' ...
                                    num2str(sum(sum(f(2,:,:))))]);
        text(3,pStill-txtOffProb,[num2str(sum(f(3,:,1))) '/' ...
                                    num2str(sum(sum(f(3,:,:))))]);
             
    % Plot probability vs. direction
    subplot(2,2,2)
        % Plot
        errorbar([1 2 3],[pToward; pSide; pAway],[iToward; iSide; iAway],'+k')
        hold on
        bar([1 2 3],[pToward; pSide; pAway],'w')
        
        % Labels & axes
        xlabel('toward                         side                         away')
        set(gca,'YTick',yTicksP)
        set(gca,'XTickLabel',' ')
        set(gca,'YLim',[yTicksP(1) yTicksP(end)]);
        hold off
        
         % Numbers
        text(1,pToward-txtOffProb,[num2str(sum(f(:,1,1))) '/' ...
                                    num2str(sum(sum(f(:,1,:))))]);
        text(2,pSide-txtOffProb,[num2str(sum(f(:,2,1))) '/' ...
                                    num2str(sum(sum(f(:,2,:))))]);
        text(3,pAway-txtOffProb,[num2str(sum(f(:,3,1))) '/' ...
                                    num2str(sum(sum(f(:,3,:))))]);
        
        
    % Plot latency vs. activity level
    subplot(2,2,3)
        % Plot
         errorbar([1 2 3]',[mean(L.a_beat); mean(L.a_glide); mean(L.a_still)],...
            [std(L.a_beat); std(L.a_glide); std(L.a_still)],'+k')
        hold on
        bar([mean(L.a_beat); mean(L.a_glide); mean(L.a_still)],'w')
        
        % Sample size
        text(1,mean(L.a_beat)+0*std(L.a_beat)-txtOffLatency,num2str(length(L.a_beat)));
        text(2,mean(L.a_glide)+0*std(L.a_glide)-txtOffLatency,num2str(length(L.a_glide)));
        text(3,mean(L.a_still)+0*std(L.a_still)-txtOffLatency,num2str(length(L.a_still)));
        
        % Labels & axes
        xlabel('beat                         glide                         still')
        set(gca,'XTickLabel',' ')
        %set(gca,'YTick',yTicksL)
        %set(gca,'YLim',[yTicksL(1) yTicksL(end)]);
        ylabel('Latency (s)')
        
        hold off
        

    % Plot latency vs. direction
    subplot(2,2,4)
        % Plot
        errorbar([1 2 3]',[mean(L.d_toward); mean(L.d_side); mean(L.d_away)],...
            [std(L.d_toward); std(L.d_side); std(L.d_away)],'+k')
        hold on
        bar([mean(L.d_toward); mean(L.d_side); mean(L.d_away)],'w')

        % Sample sizes
        text(1,mean(L.d_toward)+0*std(L.d_toward)-txtOffLatency,num2str(length(L.d_toward)));
        text(2,mean(L.d_side)+0*std(L.d_side)-txtOffLatency,num2str(length(L.d_side)));
        text(3,mean(L.d_away)+0*std(L.d_away)-txtOffLatency,num2str(length(L.d_away)));
        
        % Labels & axes
        hold off
        set(gca,'XTickLabel',' ')
        %set(gca,'YTick',yTicksL)
        %set(gca,'YLim',[yTicksL(1) yTicksL(end)]);
        xlabel('toward                         side                         away')
        
    
end

%% Run stats

if runStats
    load([mPath filesep 'ANOVA_data.mat'])
    load([mPath filesep 'probability_data.mat'])
%     % ANOVA: Test the effects on latency
%     sName = 'ltncy';
%     varnames= {['swTime_' sName] ['active_' sName] ['dir_' sName]};
%     
%     [p,table,stats] = anovan(s.latency,...
%                              {s.switchTime s.activity s.direction},'model','interaction',...
%                              'varnames',varnames);
                         
    % Chi-square: test effects of direction and activity on response
    var1 = 'activity';
    var2 = 'direction';
    var3 = 'fast start';
    
    ChiSqr3D(f,var1,var2,var3,0.05)
end


return

% Code form before:
% Calculate probability of fast start in bins 1 (swimming) and 2 (non-swimming)
prob_d1 = sum(d1)/length(d1);
prob_d2 = sum(d2)/length(d2);
prob_all = [prob_d1, prob_d2];

% Calculate the maximum likelihood (phat) and 95% confidence intervals pci for each age
[phat_1,pci_1] = binofit(sum(d1),length(d1));
[phat_2,pci_2] = binofit(sum(d2),length(d2));

figure;
bar(prob_all)
grid on;
hold on;

% Plot results
errorbar(prob_d1,phat_1',pci_1(:,1)-phat_1',pci_1(:,2)-phat_1')
errorbar(prob_d2,phat_2',pci_2(:,1)-phat_2',pci_2(:,2)-phat_2')    


