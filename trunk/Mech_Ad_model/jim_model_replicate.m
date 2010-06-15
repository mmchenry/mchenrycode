function jim_model_replicate
% This version of jim's model replicates jim's plot in the mathematica file
% "lever_afms_ndsolve_timeout.nb"

%% Parameter values set to unity

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

% Physiological cross-sectional area (m^2)
P.A = 1;

% Mass of body being moved (kg)
P.m = 1; 

% Total lever length
P.L_tot = 1; 

% Mechanical advantage
%P.MA = 1;
P.MA = [0.1 0.25 0.5 0.75 1 1.5 2]; 


%% Calculated parameters

% Input lever length
P.l_in = P.MA .* P.L_tot./(1+P.MA);

% Output lever length
P.l_out = P.L_tot ./ (1+P.MA);

% Max force produced by isometric contraction (N)
P.F_max = P.sigma_max .* P.A; 

% Hill's coefficient "b", defined for contraction from L_0
P.b = (P.a .* P.F_max .* P.L_0 .* P.Edot_max) ./ (P.A .* P.sigma_max);

% Non-dimensional mass (linear)
%P.m_star = (P.m .* P.L_0 .* P.Edot_max^2) ./ (P.A .* P.sigma_max .* P.E_max);

% Non-dimensional mass (lever)
P.m_star = (1./P.MA.^2) .* ...
           (P.m .* P.L_0 .* P.Edot_max^2) ...
           ./ (P.A .* P.sigma_max .* P.E_max);

% Non-dimension version of inverse of Hill's "a"
P.alpha = P.F_max ./ (P.a .* P.F_max);

% Curve fit parameter
P.n = exp( (3 .* P.alpha + 2) ./ (2 .* P.alpha + 2) );


%% Override default values 

P.alpha = 3;


%% Numerical simulations that vary m_star

% Time span for execution
Tspan = [0 10^8];

% Set ODE solver options
options = odeset('Events',@events);

% Initial strain 
IC(1) = 0;

% Initial strain rate;
IC(2) = 0;

% Loop through values of m_star
for i = 1:length(P.MA)
    
    % Solve the governing equation
    [t,e] = ode45(@(t,e) myode(t,e,P.m_star(i),P.alpha),Tspan,IC,options);
       
    % Store time series results
    time{i}       = t;
    ep_star{i}    = e(:,1);
    epDot_star{i} = e(:,2);
    v{i}          = e(:,2) .* (1./P.MA(i));
    lText{i}      = num2str(P.MA(i));
    
    % Stote major features
    timeMax(i) = t(end);
    vMax(i)    = max(v{i});
    
    clear e m_star
    
    disp([num2str(i) ' of ' num2str(length(P.MA))])
end

% Plot timeseries of all sims
if 1
    figure
    for i = 1:length(time)
        
        clr = rand(1,3).*.8+.1;
        lineWidth = 1.1;

        subplot(2,1,1)
        h = plot(time{i},epDot_star{i},'k-');
        set(h,'Color',clr)
        set(h,'LineWidth',lineWidth)
        hold on
        ylabel('musc. strain rate')
        
        subplot(2,1,2)
        h = plot(time{i},v{i},'k-');
        set(h,'Color',clr)
        set(h,'LineWidth',lineWidth)
        hold on
        ylabel('mass vel.');
        xlabel('time')
        
    end
    
    subplot(2,1,1)    
    legend(lText)
end 

if 1
    figure
    
    subplot(2,1,1)
    plot(log10(P.m_star),log10(timeMax),'k-o')
    xlabel('log10(m*)')
    ylabel('log10(time to max speed)')
    %title(['y = ' num2str(c(2)) ' * x' num2str(c(1)) ]) 
    axis square
    
    subplot(2,1,2)
    plot(log10(P.m_star),log10(vMax),'k-o')
    xlabel('log10(m*)')
    ylabel('log10(vMax)')
    axis square
    
    figure
    
    subplot(2,1,1)
    plot(log10(P.MA),log10(timeMax),'k-o')
    xlabel('log10(MA)')
    ylabel('log10(time to max speed)')
    %title(['y = ' num2str(c(2)) ' * x' num2str(c(1)) ]) 
    axis square
    
    subplot(2,1,2)
    plot(log10(P.MA),log10(vMax),'k-o')
    xlabel('log10(MA)')
    ylabel('log10(vMax)')
    axis square
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



