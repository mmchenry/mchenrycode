function pool_spir(rPath)
% Collects data on all sequences the the root directory that have been 
% analyzed with compare_spir. 
%
% The root directory must contain one directory for each temperature
% investigated.  he temperature diretcories must end with 'F'

%% Parameters

spirNums = 2:5;


%% Get path of data file, load data

if nargin < 1
    rPath = uigetdir(pwd,'Select root directory');
end

if 1
    
    % Get listing of temperature directories
    a1 = dir(rPath);
    
    % Initialize indices
    idx1 = 1;
    idx3 = 1;
    
    % Step through each temp
    for i = 3:length(a1)
        
        % Proceed, if directory name ends in 'F'
        if a1(i).isdir  && strcmp(a1(i).name(end),'F')
            
            % Get temp value from directory name
            T(idx1).temp = 10*str2num(a1(i).name(1)) + ...
                           10*str2num(a1(i).name(2));
            
            % Intialize variables to be stored for each spiracle
            for k = 1:length(spirNums)
                T(idx1).spir(spirNums(k)).delay = [];
                T(idx1).spir(spirNums(k)).period = [];
                T(idx1).R2(spirNums(k)).delay = [];
                T(idx1).R3(spirNums(k)).delay = [];
                T(idx1).R4(spirNums(k)).delay = [];
                T(idx1).R5(spirNums(k)).delay = [];
                T(idx1).R6(spirNums(k)).delay = [];
            end
            
            % Get subdirectories
            a2 = dir([rPath filesep a1(i).name]);
            
            % Loop through subdirectories
            for j = 3:length(a2)
                
                % Store, if directory starts with 'R'
                if a2(j).isdir && strcmp(a2(j).name(1),'R')
                    
                    % Check for 'include' file
                    a3 = dir([rPath filesep a1(i).name filesep a2(j).name ...
                        filesep 'include_in_pool.mat']);
                    
                    % Check for 'exclude' file
                    a4 = dir([rPath filesep a1(i).name filesep a2(j).name ...
                        filesep 'exclude_from_pool.mat']);
                    
                    % Proceed if comparison_data there
                    if ~isempty(a3)
                        
                        % Load 'r' structure
                        load([rPath filesep a1(i).name filesep a2(j).name ...
                            filesep 'comparison_data.mat']);
                        
                        % Load corresponding 's' structure
                        load([rPath filesep a1(i).name filesep a2(j).name '.mat'])
                        
                        % Store values in 'T' for L-R comparisons
                        % for same spiracle number
                        rNum_str = num2str(s.spir_num);
                        if length(rNum_str)==1
                            T(idx1).spir(s.spir_num).delay(end+1)  = ...
                                mean(r.delay);
                            
                            T(idx1).spir(s.spir_num).period(end+1) = ...
                                mean([r.periodM r.periodP]);
                            
                        elseif strcmp(rNum_str(1),'2')
                            sNum = str2num(rNum_str(2));
                            T(idx1).R2(sNum).delay(end+1) = ...
                                mean(r.delay);
                            
                        elseif strcmp(rNum_str(2),'2')
                            sNum = str2num(rNum_str(1));
                            T(idx1).R2(sNum).delay(end+1) = ...
                                mean(r.delay);
                            
                        elseif strcmp(rNum_str(1),'3')
                            sNum = str2num(rNum_str(2));
                            T(idx1).R3(sNum).delay(end+1) = ...
                                mean(r.delay);
                            
                        elseif strcmp(rNum_str(2),'3')
                            sNum = str2num(rNum_str(1));
                            T(idx1).R3(sNum).delay(end+1) = ...
                                mean(r.delay);
                            
                        elseif strcmp(rNum_str(1),'4')
                            sNum = str2num(rNum_str(2));
                            T(idx1).R4(sNum).delay(end+1) = ...
                                mean(r.delay);
                            
                        elseif strcmp(rNum_str(2),'4')
                            sNum = str2num(rNum_str(1));
                            T(idx1).R4(sNum).delay(end+1) = ...
                                mean(r.delay);
                        
                        elseif strcmp(rNum_str(1),'5')
                            sNum = str2num(rNum_str(2));
                            T(idx1).R5(sNum).delay(end+1) = ...
                                mean(r.delay);
                            
                        elseif strcmp(rNum_str(2),'5')
                            sNum = str2num(rNum_str(1));
                            T(idx1).R5(sNum).delay(end+1) = ...
                                mean(r.delay);
                            
                        elseif strcmp(rNum_str(1),'6')
                            sNum = str2num(rNum_str(2));
                            T(idx1).R6(sNum).delay(end+1) = ...
                                mean(r.delay);
                            
                        elseif strcmp(rNum_str(2),'6')
                            sNum = str2num(rNum_str(1));
                            T(idx1).R6(sNum).delay(end+1) = ...
                                mean(r.delay);
                            
                        end
                        clear s r
                        
                    elseif isempty(a4)
                        notRun{idx3} = ...
                            [rPath filesep a1(i).name filesep a2(j).name];
                        idx3 = idx3 + 1;   
                        
                    end
                end
            end
            
             % Iterate temperature index
            idx1 = idx1 + 1;
        end
        
       
        
    end
    
    %save([rPath filesep 'pooled_data'],'T')
    
    clear i j a1 a2 a3 idx1
    
else
    load([rPath filesep 'pooled_data']);   
end

% Report unanalyzed
disp(' ')
disp('compare_spir not yet run on the following:')
for i = 1:length(notRun)
    disp(notRun{i})
end
disp(' ')


%% Visualize results for same spiracle number

for i = 1:length(T)
   
    for j = 1:length(spirNums)
        delay_m(j)  = mean(T(i).spir(spirNums(j)).delay);
        delay_s(j)  = std(T(i).spir(spirNums(j)).delay);
        period_m(j) = mean(T(i).spir(spirNums(j)).period);
        period_s(j) = std(T(i).spir(spirNums(j)).period);
    end

    figure;
    
    subplot(2,1,1)
    errorbar(spirNums,period_m,period_s,period_s,'ko-')
    hold on
    xL = xlim;
    h = fill([xL(1) xL(2) xL(2) xL(1)],[0 0 1/15 1/15],[.5 .5 .5]);
    set(h,'EdgeColor','none')
    set(h,'FaceAlpha',0.5)
    set(gca,'XTick',spirNums)
    ylabel('Period (s)')
    xlabel('Spiracle number')
    title(['L-R Comparison at ' num2str(T(i).temp) ' deg F'])
    
    subplot(2,1,2)
    h = errorbar(spirNums,delay_m,delay_s,delay_s,'ko-');
    set(gca,'XTick',spirNums)
    hold on
    xlabel('Spiracle number')
    ylabel('Delay (s)')
    ylim([-1 1])
    xL = xlim;
    h = fill([xL(1) xL(2) xL(2) xL(1)],[-1/30 -1/30 1/30 1/30],[.5 .5 .5]);
    set(h,'EdgeColor','none')
    set(h,'FaceAlpha',0.5)
    hold off   
    
end

return

%% Visualize results for different spiracle numbers


for i = 1:length(T)
   
    for j = 1:length(spirNums)        
        R2delay_m(j)  = mean(T(i).R2(spirNums(j)).delay);
        R2delay_s(j)  =  std(T(i).R2(spirNums(j)).delay);
        R3delay_m(j)  = mean(T(i).R3(spirNums(j)).delay);
        R3delay_s(j)  =  std(T(i).R3(spirNums(j)).delay);
        R4delay_m(j)  = mean(T(i).R4(spirNums(j)).delay);
        R4delay_s(j)  =  std(T(i).R4(spirNums(j)).delay);
        R5delay_m(j)  = mean(T(i).R5(spirNums(j)).delay);
        R5delay_s(j)  =  std(T(i).R5(spirNums(j)).delay);
        R6delay_m(j)  = mean(T(i).R6(spirNums(j)).delay);
        R6delay_s(j)  =  std(T(i).R6(spirNums(j)).delay);
    end
    
    figure
   
    subplot(5,1,1)

    errorbar(spirNums,R2delay_m,R2delay_s,R2delay_s,'ko-')
    hold on
    xlabel('Spiracle number')
    ylabel('Delay (s)')
    ylim([-1 1])
    xL = xlim;
    h = fill([xL(1) xL(2) xL(2) xL(1)],[-1/15 -1/15 1/15 1/15],[.5 .5 .5]);
    set(h,'EdgeColor','none')
    set(h,'FaceAlpha',0.5)
    hold off
    title(['Delay relative to spiracle 2 at ' num2str(T(i).temp) ' deg F'])
    
    subplot(5,1,2)
    errorbar(spirNums,R3delay_m,R3delay_s,R3delay_s,'ko-')
    hold on
    xlabel('Spiracle number')
    ylabel('Delay (s)')
    ylim([-1 1])
    xL = xlim;
    h = fill([xL(1) xL(2) xL(2) xL(1)],[-1/15 -1/15 1/15 1/15],[.5 .5 .5]);
    set(h,'EdgeColor','none')
    set(h,'FaceAlpha',0.5)
    hold off
    title(['Delay relative to spiracle 3 at ' num2str(T(i).temp) ' deg F'])
    
    subplot(5,1,3)
    errorbar(spirNums,R4delay_m,R4delay_s,R4delay_s,'ko-')
    hold on
    xlabel('Spiracle number')
    ylabel('Delay (s)')
    %ylim([-1 1])
    xL = xlim;
    h = fill([xL(1) xL(2) xL(2) xL(1)],[-1/15 -1/15 1/15 1/15],[.5 .5 .5]);
    set(h,'EdgeColor','none')
    set(h,'FaceAlpha',0.5)
    hold off
    title(['Delay relative to spiracle 4 at ' num2str(T(i).temp) ' deg F'])
    
    subplot(5,1,4)
    errorbar(spirNums,R5delay_m,R5delay_s,R5delay_s,'ko-')
    hold on
    xlabel('Spiracle number')
    ylabel('Delay (s)')
    %ylim([-1 1])
    xL = xlim;
    h = fill([xL(1) xL(2) xL(2) xL(1)],[-1/15 -1/15 1/15 1/15],[.5 .5 .5]);
    set(h,'EdgeColor','none')
    set(h,'FaceAlpha',0.5)
    hold off
    title(['Delay relative to spiracle 5 at ' num2str(T(i).temp) ' deg F'])
    
    subplot(5,1,5)
    errorbar(spirNums,R6delay_m,R6delay_s,R6delay_s,'ko-')
    hold on
    xlabel('Spiracle number')
    ylabel('Delay (s)')
    %ylim([-1 1])
    xL = xlim;
    h = fill([xL(1) xL(2) xL(2) xL(1)],[-1/15 -1/15 1/15 1/15],[.5 .5 .5]);
    set(h,'EdgeColor','none')
    set(h,'FaceAlpha',0.5)
    hold off
    title(['Delay relative to spiracle 6 at ' num2str(T(i).temp) ' deg F'])
    
end