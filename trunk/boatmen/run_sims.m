function run_sims(dPath)


%% Select path & load data

if nargin<1
    dPath = uigetdir(pwd,'select Directory that contains data directories');
    if dPath==0
        return
    end
end

% Load 'd'
%global d
load([dPath filesep 'data_for_sims'])


%% Parameter values



% Inital estimate for parameter values
Cd_app0 = 7;
Cd_body0 = 1;

% Water density (kg m^-3)
rho_water = 1000;

% Viscosity of water (Pa s)
mu = 0.001;

% Number of strokes
num_strokes = 5.2;

% Which stroke to compare with the data 
stroke_num = 4;


%% Simulate each stroke in every sequence

if 1 %isempty(dir([dPath filesep 'data_for_comparison.mat']))
    
    for i = 1:length(d)
        
        k = 1;
        
        com(i).t = d(i).t;
        com(i).U = d(i).U;
        
        disp(' ')
        disp([num2str(i) ' of ' num2str(length(d)) ' sequences started ...'])
        
        for j = 1:length(d(i).pwr)
            
            % Find indices for stroke if pwr preceeds rtrn
            if (length(d(i).rtrn)>=j) && (d(i).pwr(j).idx(1) < d(i).rtrn(j).idx(1))
                idx = [d(i).pwr(j).idx d(i).rtrn(j).idx];
                rtrn_P = range(d(i).rtrn(j).t);
                run_sim = 1;
                
                % Find indices for stroke if pwr follows rtrn
            elseif (length(d(i).rtrn)>=j+1) && (d(i).pwr(j).idx(1) < d(i).rtrn(j+1).idx(1))
                idx = [d(i).pwr(j).idx d(i).rtrn(j+1).idx];
                rtrn_P = range(d(i).rtrn(j+1).t);
                run_sim = 1;
                
                % Skip if no rtrn exists for current power stroke
            else
                run_sim = 0;
                
            end
            
            if run_sim
                
                % Define model parameter values to model current stroke
                dM = d(i).pwr(j);
                dM.num_strokes = num_strokes;
                dM.stroke_num  = stroke_num;
                dM.body_len = d(i).body_len;
                dM.app_len = d(i).app_len;
                dM.body_mass = d(i).body_mass;
                dM.rtrn_P = rtrn_P;
                dM.t = d(i).t;
                dM.U = d(i).U;
                
                % Calculate initial speed
                dM.initial_speed = mean(d(i).U(idx(1:3)));
                
                % Find optimal drag coefficients
                [Cd_app,Cd_bod] = run_optimization(dM,[Cd_app0 Cd_body0],idx);
                
                
                % Store
                com(i).strk(k).idx = idx;
                com(i).strk(k).t = d(i).t(idx(1:end-1));
                com(i).strk(k).U = d(i).U(idx(1:end-1));
                com(i).strk(k).d = dM;
                com(i).strk(k).Cd_app = Cd_app;
                com(i).strk(k).Cd_bod = Cd_bod;
                
                k = k+1;
                
                clear rtrn_P
                
            end
            
            disp(['    ' num2str(j) ' of ' num2str(length(d(i).pwr)) ...
                                                ' strokes completed'])
        end
    end
    
    save([dPath filesep 'data_for_comparison'],'com')
    
    clear j i idx
    
else
    disp(' ')
    disp('Loading data_for_comparison.mat . . .')
    load([dPath filesep 'data_for_comparison'])
end


%% Visualize results

for i = 1:length(com)
    
    figure
    
    t_all = [];
    uM_all = [];
    uK_all = [];
    
    for j = 1:length(com(i).strk)
        curr = com(i).strk(j);
        r = boatmen_model(curr.d,curr.Cd_app,curr.Cd_bod);
        [tK,uM,uK] = comp_stroke(r,curr.d,curr.idx,curr.d.stroke_num);
        
        % Store
        t_all    = [t_all; tK+curr.t(1)];
        uM_all   = [uM_all; uM];
        uK_all   = [uK_all; uK];
        cd_app_tmp(j)  = curr.Cd_app;
        cd_bod_tmp(j)  = curr.Cd_bod;
        
        % Plot stroke data
        subplot(2,length(com(i).strk),j)
        h = plot(tK,uM,'r',tK,uK,'-');
        title(['Stroke ' num2str(j)])
        
        clear curr
    end
    
    body_len = com(1).strk(1).d.body_len;
    app_len  = com(1).strk(1).d.app_len;
    
    Re_bod(i)= mean(uK_all).*d(i).body_len.*rho_water ./ mu;
    Cd_app(i) = max(cd_app_tmp);
    Cd_bod(i) = max(cd_bod_tmp);
    
    %subplot(2,length(com(i).strk),j+1:j+length(com(i).strk))
    %plot(t_all,uM_all,'r',t_all,uK_all,'b')
    
end

[Re_bod,idx] = sort(Re_bod);

figure;
subplot(2,1,1)
plot((Re_bod),(Cd_app(idx)),'o-')
xlabel('log10 Re of body')
ylabel('log10 Cd of appendage')

subplot(2,1,2)
plot((Re_bod),(Cd_bod(idx)),'o-')
xlabel('log10 Re of body')
ylabel('log10 Cd of body')




return

 % To evaluate results, run simulation with optimal values
            
    
    

subplot(1,length(d(i).pwr))
            plot(tK,uK,'r',tK,uM,'-')
    
%    % Calculate mean spd from measurements
%    avg_spd = mean(d(i).U);
%     
%    % reynold number for body
%    Re_bod(i) = avg_spd.*d(i).body_len.*rho_water ./ mu;
%    [Cd_app(i),Cd_bod(i)] = run_optimization(d(i),avg_spd,[Cd_app0 Cd_body0]);
%    
%    disp(' ')
%    disp(['Finished ' num2str(i) ' of ' num2str(length(d))])
%    
% end
% 
% [Re_bod,idx] = sort(Re_bod);
% 
% figure;
% subplot(2,1,1)
% plot(log10(Re_bod),log10(Cd_app(idx)),'o-')
% xlabel('log10 Re of body')
% ylabel('log10 Cd of appendage')
% 
% subplot(2,1,2)
% plot(log10(Re_bod),log10(Cd_bod(idx)),'o-')
% xlabel('log10 Re of body')
% ylabel('log10 Cd of body')

end


function [Cd_app,Cd_bod] = run_optimization(d,X0,idx)

global dT

dT = d;

% Set options
options = optimset('OutputFcn', @outfun);

% Find optimal drag coefficient values
[X,fval] = fminsearch(@model_error,X0,options);

% Extract optimal values
Cd_app = X(1);
Cd_bod = X(2);

    function e = model_error(X)
        
        
        % Define Cd values from input vector
        Cd_app  = X(1);
        Cd_body = X(2);
        
        % Run simulation
        r = boatmen_model(d,Cd_app,Cd_body);
        
        % Evaluate model at same times as the kinematics
        [tK,uM,uK] = comp_stroke(r,d,idx,d.stroke_num);
        
        % Visualize model/data comparison 
        if 0
            figure;
            plot(tK,uK,'r',tK,uM,'-')
        end
        
        % Total error
        e = sum((uM-uK).^2);
        
    end

end

function stop = outfun(X,optimValues, state)

global dT

stop = false;

r = boatmen_model(dT,X(1),X(2));

iter = optimValues.iteration;
if 0 %isfield(d,'last_iter') && ((floor(d.last_iter/5) > floor(iter/5)))
    %figure;
    %plot(r.t,r.U,'r')
    [tK,uM,uK] = comp_stroke(r,dT,dT.idx,dT.stroke_num);
    
    plot(tK,uM,'-',tK,uK,'r-')
    title([num2str(iter) '   Cd app = ' num2str(X(1)) '   Cd body = ' num2str(X(2))])
    grid on
    pause(.2)
    %xlim([d.stroke_num.*d.spd_P*2 (d.stroke_num+1).*d.spd_P*2]-d.spd_phs)
end
d.last_iter = iter;

end

function [tK,uM,uK] = comp_stroke(r,d,idx,stroke_num)

% Extract current stroke kinematics
tK = d.t(idx(1:end-1));
uK = d.U(idx(1:end-1));

% Time vector for kinematics
tK = tK -  tK(1) + mean(diff(tK));

% Calculate index for model values after initial strokes
idxM = (r.t > (stroke_num-1).*(d.spd_P + d.rtrn_P));

tM = r.t(idxM) - r.t(find(idxM,1,'first'));
uM = r.U(idxM);

if max(tM) < max(tK)
    error('need to run simulation for more beats')
end

% Interpolate
uM = interp1(tM,uM,tK);
end


function [tK,uM,uK] = mod_compare(r,d)

% Time vector for kinematics
tK = d.t(1:end-1) -  d.t(1) + mean(diff(d.t));

% Calculate idx for measurements at start of first power stroke
idxK = tK >= d.t(1);

% Trim speed and time vectors
tK = tK(idxK) - tK(find(idxK,1,'first'));
uK = d.U(idxK);

% Calculate idx for model values after initial stroke
idxM = (r.t > (d.spd_P/2 + d.rtrn_P));

tM = r.t(idxM) - r.t(find(idxM,1,'first'));
uM = r.U(idxM);

if max(tM) < max(tK)
    error('need to run simulation for more beats')
end

% Interpolate
uM = interp1(tM,uM,tK);
end



