function d = locust_model(p,visData,visAnimate)
% Numerical solutions to locust leg kicking model
% The angular position of the lever, theta, equals zero at rest
%
% Parmeter values ('leg' = tibia + tarsus)
%   leg mass = 24e-6 (Bennet-Clark, 1975)
%   leg diameter = 1.18e-3; (Bennet-Clark, 1975)
%   leg length = 30e-3; (Bennet-Clark, 1975)
%   k spring (estimate) = 1.2e4; (Bennet-Clark, 1975)
%   Shortening length = 1.02e-3; (Bennet-Clark, 1975)
%   In-lever length = 0.76e-3 (Heitler, 1974)
%   Initial angle btwn tendon & leg = 6 deg (Heitler, 1974)

if nargin < 2
    visData    = 1;
end

if nargin < 3
    visAnimate = 0;
end


%% Default parameters

if nargin < 1
    
    % Mass (kg)
    p.m = 21e-6;
    
    % Spring stiffness (N/m)
    p.k = 2e4;
    
    % Spring stretch (m)
    p.L_s = 0.78e-3;
    
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
end


%% Calculated parameters

% Position of base of spring
p.pSpring = [p.hIn.*cos(p.theta0)-p.L_s p.hIn.*sin(p.theta0)];

% Moment of inertia (kg m^2)
%p.I = (p.m/12) * (p.leg_len^2 + p.leg_dia^2);
p.I = (p.m/3) * (p.leg_len^2);

% Moment of water inertia (kg m^2)
p.I_water = pi/3 * p.rho * (p.leg_dia/2)^2 * p.leg_len^3;

% Length of out-lever (m)
p.hOut = sqrt(p.I/p.m);

% Drag product
p.drag_prod = (pi/4) * p.Cd * p.rho * (p.leg_dia/2) * p.leg_len^4;


%% Non-dimensionalize parameters

% Constant to non-dimensionalize time
t_const = 1/sqrt(p.k/p.m);

% Specify global parameters that factor into governing equation
global I I_water L_s hIn hOut k pSpring g drag_prod m leg_len max_theta

% Non-dimensionalize parameter values
m         = p.m/p.m;
I         = p.I/p.m/(p.hIn+p.hOut)^2;
I_water   = p.I_water/p.m/(p.hIn+p.hOut)^2;
drag_prod = p.drag_prod/p.m/(p.hIn+p.hOut)^2;
k         = p.k/p.m*t_const^2;
L_s       = p.L_s/(p.hIn+p.hOut);
hIn       = p.hIn/(p.hIn+p.hOut);
hOut      = p.hOut/(p.hIn+p.hOut);
pSpring   = p.pSpring/(p.hIn+p.hOut);
g         = p.g/(p.hIn+p.hOut)*t_const^2;
leg_len   = p.leg_len/(p.hIn+p.hOut);
max_theta = p.max_theta;


%% Simulation parameters

tspan      = [0 6*10^-3]./t_const;
init_vals  = [p.theta0; 0];
options    = odeset('Events',@evnts,'RelTol',1e-7);


%% Run simulation 

[t,X] = ode45(@gov_eqn,tspan,init_vals,options);

theta  = X(:,1);
Dtheta = X(:,2);


% Clear variables
clear X m k L_r hIn hOut pSpring

% Add dimensions
Dtheta = Dtheta ./ t_const;
t      = t .* t_const;

% % Make sure sim ran for long enough
% if abs(theta(end)) > pi/100
%     error('Evaluation time not set long enough')
% end 

% Position of center of intertia of lever
pOut = [-p.hOut*cos(theta) -p.hOut*sin(theta)];

% Position of COM of lever
pCOM = [-p.leg_len/2*cos(theta) -p.leg_len/2*sin(theta)];

% Position of spring-lever intersection
pIn  = [ p.hIn*cos(theta)   p.hIn*sin(theta)];

% Length of spring
L  = pIn(:,1)-p.pSpring(1);

% Kinetic energy
E_kinetic = 0.5*(p.I+p.I_water).*Dtheta.^2;

% Elastic energy
E_potential = 0.5.*p.k.*L.^2;

% Gravitational potential energy
E_g = p.g*p.m.*(pCOM(:,2)-min(pCOM(:,2)));
tau_g = cross([pCOM zeros(length(pCOM),1)],...
        [zeros(length(pCOM),1) -p.m*p.g.*ones(length(pCOM),1) zeros(length(pCOM),1)]);
%E_g   = cumtrapz(abs(theta(1)-theta),abs(tau_g(:,3)));

% Drag torque
tau_drag = -p.drag_prod .* Dtheta .* abs(Dtheta);

% Drag energy
E_drag = cumtrapz([0; cumsum(abs(diff(theta(1)-theta)))],abs(tau_drag));

% Total energy
E_tot = E_kinetic + E_potential + E_g + E_drag;

% Re
Re = p.rho * p.leg_dia .* Dtheta .* p.leg_len ./ p.mu;

% Effective mechanical advantage
EMA = sin(theta).*(p.hIn/p.hOut);

% Velocity of center of inertia
v = Dtheta.*p.hOut;

clear t_const

disp(' ')
disp(['vmax = ' num2str(max(v))])
disp(['Prediction = ' num2str(sqrt(p.k/p.m)*p.L_s)])
disp(' ')


%% Analytical approximation

const       = sin(20/180*pi);

% theta_ana =  pi - p.pSpring(1)/p.hIn + ((p.pSpring(1)+p.hIn*(-pi+p.theta0))./p.hIn )...
%             .*cos(sqrt(const*p.k/p.m)*(p.hIn/p.hOut).*t);
%         
% Dtheta_ana = -sqrt(const*p.k/p.m) * (1/p.hOut) * (p.pSpring(1)+p.hIn*(-pi+p.theta0)) ...
%             .* sin(sqrt(const*p.k/p.m)*(p.hIn/p.hOut).*t);
%         

theta_ana =  p.L_s/p.hIn + (-p.L_s + p.hIn*p.theta0)/p.hIn ...
            .* cos(sqrt(const*p.k/p.m)*(p.hIn/p.hOut).*t);
        
Dtheta_ana = sqrt(const*p.k/p.m) .* (p.L_s - p.hIn.*p.theta0)/p.hOut ...
            .* sin(sqrt(const*p.k/p.m)*(p.hIn/p.hOut).*t);


%% Animate simulation

if visAnimate
    f = figure;
    set(f,'DoubleBuffer','on')
    
    for i = 1:length(t)
        plot([p.pSpring(1) pIn(i,1)],[p.pSpring(2) pIn(i,2)],'r',...
            [pIn(i,1) pOut(i,1)],[pIn(i,2) pOut(i,2)],'k');
        title(['t = ' num2str(t(i))])
        grid on
        xlim([-5e-3 5e-3])
        axis equal
        pause(.5)
        
    end
    
    close
end


%% Plot results

if visData
    figure;
    subplot(4,1,1)
    plot(t,theta.*(180/pi),'-',t,theta_ana.*(180/pi),'r--')
    grid on
    xlabel('Time (s)')
    ylabel('Theta (deg)')
    legend('numerical','analytical','Location','NorthWest')
    
    subplot(4,1,2)
    plot(t,Dtheta.*(180/pi/1000),'-',t,Dtheta_ana.*(180/pi/1000),'r--')
    xlabel('Time (s)')
    ylabel('dTheta/dt (deg/ms)')
    grid on
    
    subplot(4,1,3)
    plot(t,v,'-')
    xlabel('Time (s)')
    ylabel('Speed (m/s)')
    grid on
    
    subplot(4,1,4)
    plot(t,E_kinetic.*1000,'-',t,E_potential.*1000,'r',...
        t,E_g.*1000,'g',t,E_drag.*1000,'m',...
        t,E_tot.*1000,'k--')
    grid on
    ylabel('Energy (mJ)')
    xlabel('Time (s)')
    legend('kinetic','elastic','grav.','drag','total','Location','West')
    
%     figure;
%     subplot(3,1,1)
%     plot(t,Re)
%     ylabel('Re')
%     xlabel('time (s)')
%     
%     subplot(3,1,2)
%     plot(t,1./EMA)
%     ylabel('1/EMA')
%     xlabel('time (s)')
end


%% Return results

d.t = t;
d.p = p;
d.theta = theta;
d.Dtheta = Dtheta;
d.theta_ana = theta_ana;
d.Dtheta_ana = Dtheta_ana;
d.v = v;
d.EMA = EMA;

if d.theta(end) >= p.max_theta
    d.hit_max = true;
else
    d.hit_max = false;
end


function dX = gov_eqn(t,X)
global I I_water hIn hOut k pSpring g drag_prod m leg_len max_theta

% Angle of lever system
theta  = X(1);

% Rate of angular rotation of lever system
Dtheta = X(2);

% Position of intersection of spring and lever
pIn  = [ hIn*cos(theta)   hIn*sin(theta)];

% Position of center of intertia
pOut = [ -hOut*cos(theta)   -hOut*sin(theta)];

% Position of center of mass
pCOM = [-leg_len/2*cos(theta) -leg_len/2*sin(theta)];

% Current length of spring
%L  = max([pIn(1)-pSpring(1) 0]);
L  = pIn(1)-pSpring(1);

% Torque generated by spring
%tau_in = hIn*Fin.*sin(theta);
tau_spring = cross([pIn 0],[-k*L 0 0]);

% Torque from weight
tau_g = cross([pCOM 0],[0 -m*g 0]);

% Torque from drag
tau_drag = -drag_prod .* Dtheta .* abs(Dtheta);

% Define output: rate of rotation
dX(1,1) = Dtheta;

% Define output: rotational acceleration
dX(2,1) = (tau_spring(3) + tau_g(3) + tau_drag)/(I + I_water);



function [value,isterminal,direction] = evnts(t,X)
global max_theta
% Halts execution of the model

global pSpring hIn 

isterminal = 1;

%if hIn*cos(X(1))<pSpring(1)
if X(2)<=0
    value = 1;
    direction = 1;
elseif X(1) >=max_theta
    value = 0;
    direction = 0;
else
    value = 1;
    direction = 1;
end