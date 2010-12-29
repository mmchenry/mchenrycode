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
%load([dPath filesep 'data_for_sims'])
load([dPath filesep 'pooled_data'])


%% Parameter values



% Inital estimate for parameter values
Cd_app0 = 1;
Cd_body0 = 1;

% Water density (kg m^-3)
rho_water = 1000;

% Viscosity of water (Pa s)
mu = 0.001;

% Number of strokes
num_strokes = 5.2;

% Which stroke to compare with the data 
stroke_num = 4;

% Interval between plotting steps of iterations
plot_interval = 5;


%% Simulate each stroke in every sequence

if 1 %isempty(dir([dPath filesep 'data_for_comparison.mat']))
    for i = 1 %1:length(d)
        
        disp(['Starting ' num2str(i) ' of ' num2str(length(d)) ...
              '  individuals . . .'])
        
        
        dC = d(i);
        dC.num_strokes = num_strokes;
        dC.plot_interval = plot_interval;
        dC.iteration = [];
        dC.error     = [];
        
        for j = 1 %1:size(d(i).t_pwr,1)
                
                dC.t_events  = d(i).t_pwr(j,:);
                dC.U_initial = fnval(dC.sp_spd_bod,dC.t_events(1));
                
                % Find optimal drag coefficients
                [Cd_app,Cd_bod,dC] = run_optimization(dC,[Cd_app0 Cd_body0]);
                
                % Run simulation with coefficients
                r = boatmen_model(dC,Cd_app,Cd_bod);
                    
                % Plot results
                if 0
                    uK = fnval(dC.sp_spd_bod,r.t);
                    
                    figure;
                    plot(r.t,r.U,'r',r.t,uK,'-')
                end
                
                % Store
                com(i).strk(j).r = r;
                com(i).strk(j).d = dC;
                com(i).strk(j).Cd_app = Cd_app;
                com(i).strk(j).Cd_bod = Cd_bod;
                
                clear r Cd_app Cd_bod
                
                disp(['    ' num2str(j) ' of ' num2str(size(d(i).t_pwr,1)) ...
                                                ' strokes completed'])
        end
        
        clear dC
    end
    
    save([dPath filesep 'data_for_comparison'],'com')
    
    clear j i idx
    
else
    disp(' ')
    disp('Loading data_for_comparison.mat . . .')
    load([dPath filesep 'data_for_comparison'])
end

i=1;j=1;figure;plot(com(i).strk(j).d.iteration,com(i).strk(j).d.error,'o-')
return

%% Visualize results

%seqs = [1 2 4:length(com)];
seqs = [1:length(com)];

for i = 1 %seqs
    if 1
        figure
    end
    
    t_all = [];
    uM_all = [];
    uK_all = [];
    
    for j = 1 %:length(com(i).strk)
        c  = com(i).strk(j);
        uK = fnval(c.d.sp_spd_bod,c.r.t);
        
        % Store
        cd_app_tmp(j)  = c.Cd_app;
        cd_bod_tmp(j)  = c.Cd_bod;
        
        % Plot stroke data
        if 1
            subplot(2,length(com(i).strk),j)
            h = plot(c.r.t,c.r.U,'r',c.r.t,uK,'-');
            title([num2str(j) ': Cbd = ' num2str(c.Cd_bod) ...
                ' Cap= ' num2str(c.Cd_app)])
        end
        
        clear c
    end
    
    body_len = com(1).strk(1).d.body_len;
    app_len  = com(1).strk(1).d.app_len;
    U = fnval(d(i).sp_spd_bod,d(i).t);
    
    Re_bod(i) = mean(U).*d(i).body_len.*rho_water ./ mu;
    Cd_app(i) = max(cd_app_tmp);
    Cd_bod(i) = max(cd_bod_tmp);
    
    
end

return

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






 

end


function [Cd_app,Cd_bod,dT] = run_optimization(d,X0)

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
        uK = fnval(d.sp_spd_bod,r.t);
        
        % Visualize model/data comparison 
        if 0
            figure;
            plot(r.t,r.U,'r',r.t,uK,'-')
        end
        
        % Total error
        e = sum((r.U-uK).^2);
         
    end

end

function stop = outfun(X,optimValues, state)

global dT

stop = false;

r = boatmen_model(dT,X(1),X(2));

iter = optimValues.iteration;

% Evaluate model at same times as the kinematics
uK = fnval(dT.sp_spd_bod,r.t);

dT.iteration = [dT.iteration; iter];
dT.error     = [dT.error; sum((r.U-uK).^2)];

if 0 %floor(iter./dT.plot_interval) == ceil(iter./dT.plot_interval) 
    
    figure;
    plot(1000.*r.t,1000.*r.U,'r',1000.*r.t,1000.*uK,'-')
    title([num2str(iter) '   Cd app = ' num2str(X(1)) '   Cd body = ' num2str(X(2))])
    %grid on
    pause(.01)
    axis square
    ylim([0 100])
    set(gca,'YTick',1000.*[0:.02:.1])
end
end


