function force_analysis
% Analyzes force recordings for selected sequence, saves data

% Prompt for folder
if nargin < 1
    currPath = '/Volumes/data_commuter/Projects/Patek_project/force_data';
    %currPath = uigetdir(pwd,'Select folder');
end

% Get list of sequences
a = dir([currPath filesep 'force.*']);

if isempty(a)
    error('No force files!')
end


%% Load force data

% Look for data file
tmp = dir([currPath filesep 'force_pooled.mat']);

if 1 %isempty(tmp)
    
    for i = 1:length(a)
        
        % Get individual number
        iNum = str2double(a(i).name(end-7:end-5));

        % Read force text file
        if iNum>100
            % Load data from text file
            [n,Tz,Fz,Txy,Fx,Fy] = textread([currPath filesep a(i).name],...
                '%q %f %f %f %f %f','headerlines',1);
            force_mag = (Fx.^2 + Fy.^2 + Fz.^2).^0.5;
        else
            % Load data from text file
            [n,Tz,Fz,Txy,Fx,Fy] = textread([currPath filesep a(i).name],...
                '%q %f %f %q %q %q','headerlines',1);
            force_mag = Fz;
        end
        
        % Store force data
        fd(i).indiv    = iNum;
        fd(i).F        = force_mag;
        fd(i).Fz       = Fz;
        fd(i).Fx       = Fx;
        fd(i).Fy       = Fy;
        fd(i).t        = Tz;
        fd(i).filename = a(i).name;
        
        % Clear variables
        clear n Txy Fx Fy Fz Tz Txy tLen cLen iNum force
        
        % Update status
        disp(['Done ' num2str(i) ' of ' num2str(length(a))])
    end
    
    % Save data
    save([currPath filesep 'force_pooled'],'fd')
    
else
    % Load fd
    load([currPath filesep 'force_pooled.mat'])
    
end

clear tmp


%% Prompt to select duration for force analysis

figure;

for i = 1:length(fd)
    
    % Look for interval data file
    tmp = dir([currPath filesep 'intrvl_' fd(i).filename '.mat']);
    
    if isempty(tmp)
        
        disp('Select a time before impact and a time before a cavitation peak');
        disp( ' ')
        
        % Load force data
        t = fd(i).t;
        F = fd(i).F;
        
        % Display data      
        plot(t,F,'k')
        xlabel('Time')
        ylabel('Force')
        title(['Sequence ' num2str(i) ' of ' num2str(length(fd))])
        grid on
        hold on
        
        % Zoom in
        title('Zoom around duration to analyze and press return')
        zoom on
        pause
        
        % Prompt for duration to analyze
        title('Select duration to analyze')
        [x1,y1,b] = ginput(1);
        h = plot(x1,y1,'r+');
        [x2,y2,b] = ginput(1);
        delete(h)
        h = plot([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'r-');
        pause(0.5)
        hold off
   
        delete(h)
        
        t_start  = min([x1 x2]);
        t_end    = max([x1 x2]);      

        % Save data
        save([currPath filesep 'intrvl_' fd(i).filename '.mat'],'t_start','t_end')
        
        clear t F h idx t_start t_end
        
    end
end


%% Prompt to select data from initial force peak

figure;

for i = 1:length(fd)
    
    % Look for interval data file
    tmp = dir([currPath filesep 'peak_' fd(i).filename '.mat']);
    
    if isempty(tmp)
        
        % Load force data
        t = fd(i).t;
        F = fd(i).F;
        
        % Load t_start t_end
        load([currPath filesep 'intrvl_' fd(i).filename '.mat'])
        
        % Correct for baseline before strike
        F = F - mean(F(t<t_start));
        
        clear t_start t_end
        
        % Display data      
        plot(t,F,'k')
        xlabel('Time')
        ylabel('Force')
        title(['Sequence ' num2str(i) ' of ' num2str(length(fd))])
        grid on
        hold on
        
        % Zoom in
        title('Zoom into initial peak and press return')
        zoom on
        pause
        
        % Prompt for duration to analyze
        title('Select duration to analyze')
        [x1,y1,b] = ginput(1);
        h = plot(x1,y1,'r+');
        [x2,y2,b] = ginput(1);
        delete(h)
        h = plot([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'r-');
        pause(0.5)
        hold off
   
        delete(h)
        
        pk_start  = min([x1 x2]);
        pk_end    = max([x1 x2]);      

        % Save data
        save([currPath filesep 'peak_' fd(i).filename '.mat'],...
            'pk_start','pk_end')
        
        clear t F h idx pk_start pk_end
        
    end
end




%% Calculate momentum 

visSeq = 1;

if visSeq
    figure
end

for i = 1:length(fd)
    
    t = fd(i).t;
    F = fd(i).F;
    
    % Load t_start t_end
    load([currPath filesep 'intrvl_' fd(i).filename '.mat'])
    
    % Correct for baseline before strike
    F = F - mean(F(t<t_start));
    
    % Trim data after end of impact
    F = F(t<t_end);
    t = t(t<t_end);
    
    % Cumulative momentum transferred during strike
    cum_P = cumtrapz(t(t>t_start),F(t>t_start));
    
    t_P = t(t>t_start)-t_start;
    t   = t-t_start;
    
    if visSeq
        figure
        subplot(2,1,1)
        plot(t(t>0),F(t>0))
        %title(['Sequence ' num2str(i) ' of ' num2str(length(fd))])
        title(['Individual ' num2str(fd(i).indiv)])
        
        subplot(2,1,2)
        plot(t_P,cum_P)
        xlabel('time (s)')
        ylabel('cumulative momentum')
        title(['Total momentum = ' num2str(cum_P(end))])
        pause(1)
    end
    
    %Store data
    [tLen,cLen,bMass(i),k(i)] = indiv_data(fd(i).indiv);
    P_tot(i) = cum_P(end);

    clear tLen cLen F t cum_P  t_P
end

figure
% subplot(1,2,1)
% plot(bMass,P_tot,'o')
% xlabel('mass')
% ylabel('Momentum')
% axis square
% 
% subplot(1,2,2)
plot(k,P_tot,'o')
xlabel('k spring')
ylabel('Momentum')
axis square



        
function [tLen,cLen,bMass,k] = indiv_data(iNum)
% Lookup table for body size 

% Data from scaling spreadsheet
indiv_num = [3 12 120 121 122 123 125 129 132 133 137 2 5 10 16 126 128];
tot_len = [53.88 58.60 63.4 66.64 65.44 67.74 63.73 63.82 66 59 56.07 ...
           52.57 47.22 54.76 32.58 54.31 64.84];
cara_len = [14 14.52 15.457 16.877 15.59 16.24 15.89 16.427 17.3 15.31 ...
            13.93 13.45 11.39 14.22 7.83 14.13 15.95];
b_mass = [3.58 4.2 5.47 5.36 5.42 4.78 3.58 4.2 6.58 6.497 5.82 6.147 ...
          4.78 3.89 3.13 2.08 3.59 0.7 3.69 5.52];
    
k_spring = [48.9023 50.4358 27.32212 23.27122 38.2438 30.59903 38.65258 ...
            32.3173 33.64915 72.30235 23.12187 37.15923 40.5199 49.32086 ...
            nan nan nan nan nan nan];
        
% Find match for individual
indx = find(indiv_num==iNum);

% Return matching values
if isempty(indx)
    tLen = nan;
    cLen = nan;
    bMass = nan;
    k_spring = nan;
else
    tLen = tot_len(indx);
    cLen = cara_len(indx);
    bMass = b_mass(indx);
    k = k_spring(indx);
end




