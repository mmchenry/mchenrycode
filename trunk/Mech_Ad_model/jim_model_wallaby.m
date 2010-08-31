function jim_model



%% Define paths

dPath = '/Volumes/Docs/Projects/Velocity advantage';


%% Loaded parameters

%P = unity_params;
P = sim_params;


%% Override default values (pre calculated params)

% Elevate mass to ground reaction force for hopping
%P.m = P.m*6; 

% Vary MA
%P.MA = P.MA .* 10.^linspace(-2,2,10); 

%% Calculated parameters

% Hill's coefficient "b", defined for contraction from L_0
P.b = (P.a .* P.F_max .* P.L_0 .* P.Edot_max) ./ (P.A .* P.sigma_max);

% Non-dimension version of inverse of Hill's "a"
P.alpha = P.F_max ./ (P.a .* P.F_max);

% Curve fit parameter
P.n = exp( (3 .* P.alpha + 2) ./ (2 .* P.alpha + 2) );


%% Override default values (post calculated params)

%P.alpha = 3;


%% Numerical simulations that vary MA and mass


% Look for data file
a = dir([dPath filesep 'simData.mat']);

% Execute if data file not present
if isempty(a)
    
    % Time span for execution
    Tspan = [0 10^8];
    
    % Set ODE solver options
    options = odeset('Events',@events);
    
    % Initial strain
    IC(1) = 0;
    
    % Initial strain rate;
    IC(2) = 0;
    
    % Define MA and mass values to loop through
    MAVals = 10.^linspace(-3,1,10);
    massVals = 10.^linspace(-3,3,5);
    
    % Loop through mass values
    for i = 1:length(massVals)
        
        % Start timer
        tic
        
        % Define current mass
        mass = massVals(i);
        
        % Non-dimensional mass (for muscle)
        m0_star = (mass .* P.L_0 .* P.Edot_max^2) ...
            ./ (P.A .* P.sigma_max .* P.E_max);
        
        
        % Store
        d(i).mass          = mass;
        d(i).m0_star       = m0_star;
        d(i).alpha         = P.alpha;
        d(i).paramVals     = P;
        
        % Loop through values of MA
        for j = 1:length(MAVals)
            
            % Define current MA
            MA = MAVals(j);
            
            % Non-dimensional mass (for lever system)
            m_star = (1./MA.^2) .* m0_star;
            
            % Solve the governing equation
            [t,e] = ode45(@(t,e) myode(t,e,m_star,P.alpha),Tspan,IC,options);
            
            % Store simulation results
            d(i).sim(j).MA         = MA;
            d(i).sim(j).m_star     = m_star;
            d(i).sim(j).time       = t;
            d(i).sim(j).ep_star    = e(:,1);
            d(i).sim(j).epDot_star = e(:,2);
            d(i).sim(j).v          = e(:,2) .* (1./MA);
            d(i).sim(j).lText      = num2str(MA);
            
            % Update status
            disp(['   Done ' num2str(j) ' of ' num2str(length(MAVals)) ' sims'])
            
            % Clear variables
            clear e m_star t MA
        end
        
        % Update status
        T = toc/60;
        timeLeft = T*(length(massVals)-i);
        disp(['Done ' num2str(i) ' of ' num2str(length(massVals)) ...
            ' rounds (' num2str(timeLeft) ' min remaining)'])
        
        % Clear variables
        clear mass j T timeLeft
    end
    
    % Save data
    disp(' ')
    disp(['Saving to ' dPath filesep 'simData.mat'])
    save([dPath filesep 'simData.mat'],'d')
    
end


%% Calculate performance metrics from sims

% Load data
load([dPath filesep 'simData.mat']);

MAactual = d(1).paramVals.MA;

% Loop through mass data
for i = 1:length(d)
    
    % Store m0_star value
    m0_vals(i,1) = d(i).m0_star; 
    
    lText{i} = num2str(d(i).m0_star);
    
    % Loop through sims
    for j = 1:length(d(i).sim)
        
        % Store MA   
        MA(j,i) = d(i).sim(j).MA;
        
        accel = diff(d(i).sim(j).v)./diff(d(i).sim(j).time);
        
        % Store kinematic parameters
        accelMax(j,i) = max(accel);
        velMax(j,i)   = max(d(i).sim(j).v);
        t_velMax(j,i) = max(d(i).sim(j).time);
        
        clear accel
    end
end

clear d


%% Plot performance


figure

subplot(3,1,1)
loglog(MA(:,1),accelMax)
hold on
loglog([MAactual MAactual],ylim,'k-')
xlabel('MA')
ylabel('max accel')
legend(lText,'Location','NorthWest')

subplot(3,1,2)
loglog(MA(:,1),velMax)
hold on
loglog([MAactual MAactual],ylim,'k-')
xlabel('MA')
ylabel('max vel')

subplot(3,1,3)
loglog(MA(:,1),t_velMax)
hold on
loglog([MAactual MAactual],ylim,'k-')
xlabel('MA')
ylabel('time to max vel')


%% % Analytical prediction

% Max strain rate (muscle)
Edot_star = (1 - exp( -4 ./ (3 .* P.m_star * (1+P.alpha)))).^(1/P.n);

% Max velocity of mass
vAna = (1 ./ P.MA) .* ...
       (1-exp(-4./(3.*P.m_star.*(1+P.alpha)))).^(1/P.n);



return


 % Store performance parameters
        timeMax(j) = t(end);
        vMax(j)    = max(v{j});


% Plot timeseries of all sims
if 1
    figure
    for i = 1:length(time)
        
        clr = rand(1,3).*.8+.1;
        lineWidth = 1.1;
        
        subplot(3,1,1)
        h = plot(time{i},1-ep_star{i},'k-');
        set(h,'Color',clr)
        set(h,'LineWidth',lineWidth)
        ylabel('1 - musc. strain')
        hold on
        
        subplot(3,1,2)
        h = plot(time{i},epDot_star{i},'k-');
        set(h,'Color',clr)
        set(h,'LineWidth',lineWidth)
        hold on
        ylabel('musc. strain rate')
        
        subplot(3,1,3)
        h = plot(time{i},v{i},'k-');
        set(h,'Color',clr)
        set(h,'LineWidth',lineWidth)
        hold on
        ylabel('mass vel.');
        xlabel('time')
        
    end
    
    subplot(3,1,1)    
    legend(lText)
end 
    
% Plot wrt m*
if 1
    figure
        
    subplot(2,1,1)
    plot(log10(P.m_star),log10(timeMax),'k-o')
    hold on;
    plot(log10(P.m_star_real).*[1 1],ylim,'r-')
    xlabel('log10(m*)')
    ylabel('log10(time to max speed)')
    %title(['y = ' num2str(c(2)) ' * x' num2str(c(1)) ]) 
    axis square
    
    subplot(2,1,2)
    plot(log10(P.m_star),log10(vMax),'k-o',...
        log10(P.m_star),log10(vAna),'g--')
    hold on;
	plot(log10(P.m_star_real).*[1 1],ylim,'r-')
    xlabel('log10(m*)')
    ylabel('log10(vMax)')
    axis square
    
end

% Plot wrt MA
if 1
    figure
    
    subplot(2,1,1)
    plot(log10(P.MA),log10(timeMax),'k-o')
    hold on
    plot([1 1].*log10(P.MA_real),ylim,'r-')
    hold off
    xlabel('log10(MA)')
    ylabel('log10(time to max speed)')
    %title(['y = ' num2str(c(2)) ' * x' num2str(c(1)) ]) 
    axis square
    
    subplot(2,1,2)
    plot(log10(P.MA),log10(vMax),'k-o')
    hold on
    plot([1 1].*log10(P.MA_real),ylim,'r-')
    plot(log10(P.MA),log10(vAna),'g--')
    hold off
    xlabel('log10(MA)')
    ylabel('log10(vMax)')
    axis square
end


%% Analytical predictions

% Max strain rate (muscle)
Edot_star = (1 - exp( -4 ./ (3 .* P.m_star * (1+P.alpha)))).^(1/P.n);

% Max velocity of out lever
vAna = P.m_star .* (1 ./ P.MA) .* ...
    (1 - exp( (-4 .* P.MA.^2) ./ ...
    (3 .* P.m_star .* (1 + P.alpha)))).^(1/P.n);

% Plot data
if 0
    figure;
    subplot(2,1,2)
    plot(P.MA,v,'g-')
    xlabel('EMA')
    ylabel('v')
end




function eDot = myode(t,e,mStar,alpha)
% Governing equation of the model 
% e - non-dimensional strain rate (epsilon star in jim's notes)

% First derivative of strain (strain rate)
eDot(1,1) = e(2);

% Second derivative of strain
eDot(2,1) = (1 ./ mStar) .*  (1 - e(1).^2) .* ...
            ( (1 - e(2)) ./ (1 + alpha.*e(2)));
   
        
function [value,isterminal,direction] = events(t,e)
% Stops ode45 when negative or zero strain is reached

value      = max([1-e(1) 0]);
isterminal = 1;
direction  = 0;


function P = sim_params
% Sources for verts in general and wallaby's in particular
% A83 - Alexander (1983) Animal Mechanics
% W85 - Woledge, Curtin, & Homsher (1985) Energetic Aspects of Muscle
% Contraction
% BB95 - BIEWENER and Baudinette (1995) JEB 198:1829-1841
% B05 - BIEWENER. (2005) JEB 208:1665-1676
% BKB98 - Biewener et al (1998) 201:1681-1694
% M08 - McGowan et al. (2008)  J Anat. 212:153-163

% -------------------------------------------------------------------------
% The following parameters are generic vertebrate muscle parameters, 
% from W85

% Max stress for vert. muscle (at L_0, isometric contraction) (Pa)
P.sigma_max = 300e3; 

% Max strain, where muscle produces no force (isometric contraction)
% (dimensionless)
P.E_max = 0.4; %W85 (p. 40)

% Strain rate at which the muscle produces no force (isotonic contraction)
% (1/s)
P.Edot_max = 2; %W85 (p. 62)


% -------------------------------------------------------------------------
% In the Wallaby, Both the Plantaris and Gastroc produce significant force,
% as determined by BKB98 & BB95.
% We therefore used the sum of physiological cross-sectional area for both
% from BB95.

% Physiological cross-sectional area (m^2)
P.A = 23.57e-4; 

% This corresponds to a wallaby with the following avg. body mass, which we
% take as the mass that is accelerated (kg) (BB95)
P.m = 4.8;


% -------------------------------------------------------------------------
% The Plantaris and Gastroc have similar fiber lengths and angles. BB95
% reported mean fiber lengths among both of 16.9 mm and angles of 33.  We
% assume that these resting lengths may approximate that at which max force
% is generated.  The component of this force along the tendon is found by
% multiplying by cos(fiber angle)

% Muscle fiber length at which max force is produced (isometric contractions) (m)
P.L_0 = 10*cos(33*pi/180) .* 16.9e-3; 


% -------------------------------------------------------------------------
% The Hill coefficient "a" is approximately 0.25*F_max (for fast striated
% muscle, A83 (p. 107), W85)

% Max force produced by isometric contraction (N)
P.F_max = P.sigma_max .* P.A; 

P.a = 0.25 .* P.F_max;


% -------------------------------------------------------------------------
% Skeletal morphometrics were obtained from the scaling analysis of 
% wallabies by M08.

% Total lever length, taken as the length of the metatarsals (m)
P.L_tot = (39.88 .* P.m.^0.37) ./1000; 

% The mechanical advantage was calculated from the weighted calculation of
% moment arm, divided by the difference of metatarsal length and the moment
% arm
r_in = (12.17 .* P.m.^0.43) / 1000; % (m)
P.MA = r_in / (P.L_tot - r_in);     % (dimensionless)



function P = unity_params

%% Generic parameter values

% Max stress (at L_0, isometric contraction) (Pa)
P.sigma_max = 1; 

% Muscle length at which max force is produced (isometric contractions) (m)
P.L_0 = 1;

% Hill's coefficient "a", divided by F_max, defined for contraction from L_0
P.a = 1; 

% Max strain, where muscle produces no force (isometric contraction)
P.E_max = 1; 

% Strain rate at which the muscle produces no force (isotonic contraction)
P.Edot_max = 1; 


%% Particular parameter values
% Specific to hopping wallaby

% Physiological cross-sectional area (m^2)
P.A = 1;

% Mass of body being moved (kg)
P.m = 1; 

% Total lever length
P.L_tot = 1; 

% Effective mechanical advantage
P.MA = 1; 

% Input lever length (not needed for simulations)
P.l_in = MA .* P.L_tot./(1+MA);

% Output lever length (not needed for simulations)
P.l_out = P.L_tot ./ (1+MA);

% Max force produced by isometric contraction (N)
P.F_max = P.sigma_max .* P.A; 

