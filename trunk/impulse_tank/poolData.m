function poolData(mPath)

% mPath    - root path that holds the sequence directories (optional)
% 

%% Options for what to run

runCollect  = 1;
runPlotting = 1;
runStats    = 0;


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
    
    % Initiate p
    p.d_toward = [];
    p.d_side   = [];
    p.d_away   = [];
    p.a_beat   = [];
    p.a_glide  = [];
    p.a_still   = [];
    
    % Initiate L
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
                    
                    if (theta > 0) && (theta <= 60)
                        direction = 'toward';
                    elseif (theta > 60) && (theta <= 120)
                        direction = 'side  ';
                    elseif (theta > 120) && (theta <= 180)
                        direction = 'away  ';
                    end
                    
                    % Calculate latency
                    [amp,daq,vid] = synchronizeData(...
                                        [mPath filesep batch{i}.filename]);
                    preTrig_frms  = vid.dur - vid.postTrg;
                    
                    % Define start frame
                    start_frm = c(j).fast_start_frame_num;
                    
                    % If fast start
                    if start_frm~=0
                        
                        % Latency is the video time for the frame at which
                        % motion was visible - video time of pretrigger period
                        % - the delay in the satrt of the motor motion
                        latency = vid.t(start_frm+preTrig_frms) - amp.t_mStart;
                        
                        %[vid.t(preTrig_frms) start_frm latency]
                        
                        if latency < 0
                            warning(['start_frm = ' num2str(start_frm)...
                                ' in fish ' num2str(j) ' ' ...
                                batch{i}.filename ...
                                ' is before stimulus '...
                                '-- not recorded']);
                        else
                            
                            % Store data for latency ANOVA
                            s.activity(k,:)   = c(j).swimming;
                            s.switchTime(k,1) = c(j).switchTime;
                            s.direction(k,:)  = direction;
                            s.latency(k,1)    = latency;
                            
                            % Store data for probablilty & latency analysis
                            if strcmp(direction,'toward')
                                p.d_toward = [p.d_toward; 1];
                                L.d_toward = [L.d_toward; latency];
                            elseif strcmp(direction, 'side  ')
                                p.d_side = [p.d_side; 1];
                                L.d_side = [L.d_side; latency];
                            elseif strcmp(direction, 'away  ')
                                p.d_away = [p.d_away; 1];
                                L.d_away = [L.d_away; latency];
                            end
                            
                            if strcmp(c(j).swimming,'beats')
                                p.a_beat = [p.a_beat; 1];
                                L.a_beat = [L.a_beat; latency];
                            elseif strcmp(c(j).swimming,'glide')
                                p.a_glide = [p.a_glide; 1];
                                L.a_glide = [L.a_glide; latency];
                            elseif strcmp(c(j).swimming,'still')
                                p.a_still = [p.a_still; 1];
                                L.a_still = [L.a_still; latency];
                            end
                            
                            % Update k
                            k = k + 1;
                            
                        end
                        
                        % If no fast start
                    else
                        
                        % Store data for probablilty analysis
                        if strcmp(direction,'toward')
                            p.d_toward = [p.d_toward; 0];
                        elseif strcmp(direction,'side')
                            p.d_side = [p.d_side; 0];
                        elseif strcmp(direction, 'away  ')
                            p.d_away = [p.d_away; 0];
                        end
                        
                        if strcmp(c(j).swimming,'beats')
                            p.a_beat = [p.a_beat; 0];
                        elseif strcmp(c(j).swimming,'glide')
                            p.a_glide = [p.a_glide; 0];
                        elseif strcmp(c(j).swimming,'still')
                            p.a_still = [p.a_still; 0];
                        end
                        
                    end
                    
                    clear amp daq vid
                    clear direction theta latency
                    
                    %TO DO: figure out stats,
                    %        make graphs
                end
            end            
        end
    end
    
    save([mPath filesep 'ANOVA_data.mat'],'s')
    save([mPath filesep 'probability_data.mat'],'p')
    save([mPath filesep 'latency_data.mat'],'L')
    
end


%% Visualize results

if runPlotting
    
    % Load data
    load([mPath filesep 'latency_data.mat'])
    load([mPath filesep 'probability_data.mat'])
    load([mPath filesep 'ANOVA_data.mat'])
    
    % Parameters for plotting
    yTicksL = (floor(min(1000.*s.latency)):10:ceil(max(1000.*s.latency)))./1000;
    yTicksP = 0:.5:1;
    txtOffLatency = .004;
    txtOffProb = 0.05;
    
    % Calculate probability stats
    %pToward            = sum(p.d_toward)/length(p.d_toward);
    [pToward,iToward] = binofit(sum(p.d_toward),length(p.d_toward));
    iToward = iToward(2) - pToward;
    
    [pAway,iAway] = binofit(sum(p.d_away),length(p.d_away));
    iAway = iAway(2) - pAway;
    
    [pSide,iSide] = binofit(sum(p.d_side),length(p.d_side));
    iSide = iSide(2) - pSide;
    
    [pBeat,iBeat] = binofit(sum(p.a_beat),length(p.a_beat));
    iBeat = iBeat(2) - pBeat;
    
    [pGlide,iGlide] = binofit(sum(p.a_glide),length(p.a_glide));
    iGlide = iGlide(2) - pGlide;
    
    [pStill,iStill] = binofit(sum(p.a_still),length(p.a_still));
    iStill = iStill(2) - pStill;
    
    
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
        text(1,pBeat-txtOffProb,[num2str(sum(p.a_beat)) '/' num2str(length(p.a_beat))]);
        text(2,pGlide-txtOffProb,[num2str(sum(p.a_glide)) '/' num2str(length(p.a_glide))]);
        text(3,pStill-txtOffProb,[num2str(sum(p.a_still)) '/' num2str(length(p.a_still))]);
        
        
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
        text(1,pToward-txtOffProb,[num2str(sum(p.d_toward)) '/' num2str(length(p.d_toward))]);
        text(2,pSide-txtOffProb,[num2str(sum(p.d_side)) '/' num2str(length(p.d_side))]);
        text(3,pAway-txtOffProb,[num2str(sum(p.d_away)) '/' num2str(length(p.d_away))]);
        
        
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
        set(gca,'YTick',yTicksL)
        set(gca,'YLim',[yTicksL(1) yTicksL(end)]);
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
        set(gca,'YTick',yTicksL)
        set(gca,'YLim',[yTicksL(1) yTicksL(end)]);
        xlabel('toward                         side                         away')
        
    
end

%% Run stats

if runStats
    
    % Test the effects on latency
    sName = 'ltncy';
    varnames= {['swTime_' sName] ['active_' sName] ['dir_' sName]};
    
    [p,table,stats] = anovan(s.latency,...
                             {s.switchTime s.activity s.direction},'model','interaction',...
                             'varnames',varnames);
    
    
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


