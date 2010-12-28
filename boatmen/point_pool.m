function point_pool(dPath)


%TODO: average among power strokes to find model parameter values

% Compare measured and simulated velocities

vis_seq = 0;


%% Select path

if nargin<1
    dPath = uigetdir(pwd,'select Directory that contains data directories');
    if dPath==0
        return
    end
end


%% Pool together all sequences that have been completely analyzed

if isempty(dir([dPath filesep 'pooled_data.mat']))
    
    a = dir(dPath);
    
    idx = 1;
    
    for i = 3:length(a)
        if a(i).isdir
            
            % Check for data files
            % Look for "seq"
            tmp1 = dir([dPath filesep a(i).name filesep 'seq_info.mat']);
            
            if isempty(tmp1)
                tmp1 = dir([dPath filesep a(i).name filesep 'seq_info2.mat']);
            end
            
            % Look for "boat_coord"
            tmp2 = dir([dPath filesep a(i).name filesep 'boat_coords.mat']);
            
            % Look for "calconst"
            tmp3 = dir([dPath filesep a(i).name filesep 'cal_const.mat']);
            
            % Look for "body_length"
            tmp4 = dir([dPath filesep a(i).name filesep 'body_mass.mat']);
            
            % Issue warnings & don't run analysis, if data files not present
            if isempty(tmp1)
                warning([dPath filesep a(i).name filesep ...
                          'seq_info.mat does not exist'])
                
            elseif isempty(tmp2)
                warning([dPath filesep a(i).name filesep ...
                           'boat_coords.mat does not exist'])
         
           elseif isempty(tmp3)
                warning([dPath filesep a(i).name filesep ...
                           'cal_const.mat does not exist'])
                       
           elseif isempty(tmp4)
                warning([dPath filesep a(i).name filesep ...
                           'body_mass.mat does not exist'])
                
            % Otherwise . . .
            else
                
                % Check for complete dataset
                load([dPath filesep a(i).name filesep 'boat_coords.mat'])
                % Index of frames for which there are data for both the wrist and tip
                nonan = ~(isnan(b.xTip)  | isnan(b.xWrist) | isnan(b.xNose) | isnan(b.xTail));
                
                % Chech that frames of data exist
                if sum(nonan)==0
                    
                    warning(['Incomplete dataset in ' dPath filesep ...
                        a(i).name filesep 'boat_coords.mat'])
                    
                else
                    
                    % Get the 'd' structure
                    d(idx).seq  = a(i).name;
                    tmp = point_analyze([dPath filesep a(i).name]);
                    
                    d(idx).t           = tmp.t';
                    d(idx).frames      = tmp.frames';
                    d(idx).pwr         = tmp.pwr;
                    d(idx).rtrn        = tmp.rtrn;
                    d(idx).body_len    = tmp.body_len;
                    d(idx).app_len     = tmp.app_len;
                    d(idx).cntr_G      = tmp.cntr_G;
                    
                    % Get body mass
                    load([dPath filesep a(i).name filesep 'body_mass.mat']);
                    
                    d(idx).body_mass = bMass;
                    
                    clear bMass tmp
                    
                    
                    % Visualize results from each sequence
                    if vis_seq
                        figure
                        
                        % Step though each power stroke
                        for j = 1:length(d(idx).pwr)
                            
                            % Speed and corresponding time 
                            t = d(idx).pwr(j).t(2:end)-d(idx).pwr(j).t(1);
                            spd = d(idx).pwr(j).spd;
                            
                            %data_filtered = butter_filt(data,sample_rate,cut_freq,type) 
                            
                            % Speed (from function)
                            spd_f = spd_func(t,d(idx).pwr(j).spd_A,...
                                d(idx).pwr(j).spd_P,d(idx).pwr(j).spd_phs,...
                                d(idx).pwr(j).spd_0);
                            
                            % Angle and corresponding time
                            tA = d(idx).pwr(j).t-d(idx).pwr(j).t(1);
                            ang = d(idx).pwr(j).angl;
                            
                            % Angle (from function)
                            ang_f = ang_func(tA,d(idx).pwr(j).ang_amp,...
                                                 d(idx).pwr(j).ang_P,...
                                                 d(idx).pwr(j).ang_start);
                            
                            % Plot                  
                            subplot(2,length(d(idx).pwr),j)
                            h = plot(t,spd,[clrs(j) 'o'],t,spd_f,...
                                                          [clrs(j) '-']);
                            set(h(1),'MarkerfaceColor',clrs(j))
                            hold on
                            title([a(i).name ' stroke ' num2str(j)])
                            
                            subplot(2,length(d(idx).pwr),j+length(d(idx).pwr))
                            h = plot(tA,ang,[clrs(j) 'o'],...
                                     tA,ang_f,[clrs(j) '-']);
                            set(h(1),'MarkerfaceColor',clrs(j))
                            hold on
                            
                        end
                        
                        subplot(2,length(d(idx).pwr),1)
                        ylabel('spd')
                        
                        subplot(2,length(d(idx).pwr),1)
                        ylabel('angl')
                        
                    end
                    
                    idx = idx + 1;
                end
            end
        end
    end
    
    % Save data
    save([dPath filesep 'pooled_data'],'d')
    
else
    
    % Loading pooled data
    disp('Loading pooled_data . . .')
    load([dPath filesep 'pooled_data'])
    disp(' ')
    
end


%% Find model parameters for each sequence


for i = 1:length(d)
    
    % Step through power strokes
    for j = 1:length(d(i).pwr)
       spd_A(j)     = d(i).pwr(j).spd_A;
       spd_phs(j)   = d(i).pwr(j).spd_phs;
       spd_P(j)     = d(i).pwr(j).spd_P;
       
       ang_P(j)     = d(i).pwr(j).ang_P;
       ang_amp(j)   = d(i).pwr(j).ang_amp;
       ang_start(j) = d(i).pwr(j).ang_start;
       ang_start(j) = d(i).pwr(j).ang_start;
    end
    
    % Step through return strokes
    for j = 1:length(d(i).rtrn)
        periods(j) = d(i).rtrn(j).period;
    end
    
    % Store average speed variables for power stroke
    d(i).spd_A     = mean(spd_A);
    d(i).spd_phs   = mean(spd_phs);
    d(i).spd_P     = mean(spd_P);
    d(i).ang_P     = mean(ang_P);
    d(i).initial_speed = 0; 
    %TODO: calculate initial speed for each sequence
    
    % Store average angular variables for power stroke
    d(i).ang_amp   = mean(ang_amp);
    d(i).ang_start = mean(ang_start);
    
    % Store average return stroke period
    d(i).rtrn_P    = mean(periods);
    
    d(i).num_strokes = 3;
    
    clear spd_A spd_phs spd_P ang_P ang_amp ang_start
    
    % Calculate COM speed
    d(i).U = sqrt(diff(d(i).cntr_G(:,1)).^2 + diff(d(i).cntr_G(:,1)).^2)./diff(d(i).t);
    
    
    % Plot simulation
    if 0
        
        r(i) = boatmen_model(d(i));
        
        figure;
        subplot(3,1,1)
        [ax,h1,h2] = plotyy(r(i).t.*1000,r(i).gamma,r(i).t.*1000,r(i).v_n.*1000);
        ylabel(ax(1),'gamma (rad)')
        ylabel(ax(2),'v_n (mm/s)')
        grid on
        
        subplot(3,1,2)
        plot(r(i).t.*1000,r(i).thrust.*10^6,'g',r(i).t.*1000,r(i).drag.*10^6,'-')
        legend('thrust','drag')
        ylabel('Force (micro N)')
        grid on
        
        subplot(3,1,3)
        plot(r(i).t.*1000,r(i).U.*1000,'-',...
             (d(i).t(2:end)-d(i).t(1)).*1000,d(i).U.*1000,'r-')
    end  
    
end

% Save data
disp('')
disp('Saving data_for_sims . . .')
save([dPath filesep 'data_for_sims'],'d')




function y = spd_func(t,A,P,phs,spd_0)
% Function that defines the speed of the power stroke

y = spd_0 + A.*sin(pi.*(t-phs)./(1.5*P)).^2;


function y = ang_func(t,A,P,ang_start)
% Defines angle of power stroke

y = ang_start + A.*sin(pi.*(t)./(2.*(1.1*P))).^2;


function c = clrs(n)

cs{1} = 'r';
cs{2} = 'b';
cs{3} = 'g';
cs{4} = 'm';
cd{5} = 'k';

c = cs{n};


function data_filtered = butter_filt(data,sample_rate,cut_freq,type) 
% High-pass or low-pass butterworth filter

% All frequency values are in Hz.

% Nyquist freq.
Nqst = sample_rate/2;   

% Calculate stopband frequency
if strcmp(type,'high')
    stop_freq = max([(cut_freq - Nqst/10) .01]);  

elseif strcmp(type,'low')
    stop_freq = min([(cut_freq + Nqst/10) (Nqst-.01)]); 
 
end

% Stopband Attenuation (dB)
Astop = 10;

% Passband Ripple (dB)
Apass = 1;   

% Normalise the cutoff freq. to the Nyquist freq
Wp = cut_freq/Nqst;

% Normalise the stoppass freq. to the Nyquist freq
Ws = stop_freq/Nqst;

% Check cutoff
if (Wp > 1) || (Ws > 1)
    error('Cutoff or bandpass frequencies cannot exceed the Nyquist freq')
end

% Calculate the order from the parameters using BUTTORD.
[N,Fc] = buttord(Wp, Ws, Apass, Astop);    
    
% Calculate the zpk values using the BUTTER function.
[B A] = butter(N, Fc, type);

% Plot frequency reponse
%freqz(B,A,512,sample_rate); 

% Filter the data
data_filtered   = filtfilt(B,A,data); 

