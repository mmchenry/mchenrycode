function run_locust_sims

%TODO: Fix scaling

%% Default parameters

% Mass (kg)
p.m = 1*24e-6;

% Spring stiffness (N/m)
p.k = 1e4;

% Spring stretch (m)
p.L_s = 1e-3;

% Length of in-lever (m)
p.hIn = 0.76e-3;

% Diameter of leg (m)
p.leg_dia = 1.18e-3;

% Length of leg (m)
p.leg_len = 30e-3;

% Initial angle of lever system (rad)
p.theta0 = 6/180*pi;

% Acceleration of gravity (m/s^2)
p.g = 9.81;

% Density of air (kg m^-3)
p.rho = 1.2;

% Dynamic viscosity of air (Pa s)
p.mu = 1.845e-5;

% Drag coefficient for a cylinder (dimensionless)
p.Cd = 1.2;

% Max excursion of theta (rad)
p.max_theta = 200/180*pi;

%% Calculated parameters

% Position of base of spring
p.pSpring = [p.hIn.*cos(p.theta0)-p.L_s p.hIn.*sin(p.theta0)];

% Moment of inertia (kg m^2)
p.I = (p.m/12) * (p.leg_len^2 + p.leg_dia^2);

% Moment of water inertia (kg m^2)
p.I_water = pi/3 * p.rho * (p.leg_dia/2)^2 * p.leg_len^3;

% Length of out-lever (m)
p.hOut = sqrt(p.I/p.m);

% Drag product
p.drag_prod = (pi/4) * p.Cd * p.rho * (p.leg_dia/2) * p.leg_len^4;



%% Vary in-lever length

% Set values for in-levers
in_levers = p.hIn .* 10.^[-1:.2:1];

% Remove fluid forces and gravity
p.g = 0;
p.rho = 0;
    

for i = 1:length(in_levers)
    
    % Set current in-lever
    p.hIn = in_levers(i);
    
    d = locust_model(p,0,0);
    
    max_v(i) = max(d.v);
    max_t(i) = d.t(end-1);
    EMA(i) = max(d.EMA);
    hit_max(i) = d.hit_max;
    
    clear d
    
end

% Scaling of period
c = polyfit(log10(in_levers(hit_max)./p.hOut),log10(max_t(hit_max)),1);


figure
subplot(1,2,1)
plot(in_levers(hit_max)/p.hOut,max_v(hit_max),'o-',...
    in_levers(~hit_max)/p.hOut,max_v(~hit_max),'ro-')
xlabel('MA')
ylabel('max vel')

subplot(1,2,2)
plot(in_levers(hit_max)/p.hOut,max_t(hit_max),'o-',...
     in_levers(~hit_max)/p.hOut,max_t(~hit_max),'ro-')
xlabel('MA')
ylabel('time to max vel')

figure
plot(log10(in_levers(hit_max)./p.hOut),log10(max_t(hit_max)),'b-o',...
     log10(in_levers(hit_max)./p.hOut),...
     log10(c(2).*(in_levers(hit_max)./p.hOut).^c(1)),'g--')
title(['fit = ' num2str(c(1)) ])
%axis equal