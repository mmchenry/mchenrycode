function [d,p] = locust_model(p,visData,visAnimate)
% Numerical solutions to locust leg kicking model
% The angular position of the lever, theta, equals zero at rest


if nargin < 2
    visData    = 1;
end

if nargin < 3
    visAnimate = 0;
end


%% Default parameters

if nargin < 1
    
    % Use defult parameter values, if none given
    p = default_param;

end


%% Calculated parameters

% Position of base of spring
%p.pSpring = [p.hIn.*cos(p.theta0)-p.L_s p.hIn.*sin(p.theta0)];
%p.pSpring = [0 p.hIn.*sin(p.theta0)];


p.xJoint = p.L_s-p.hIn.*cos(p.theta0);

% Moment of inertia (kg m^2)
%p.I = (p.m/12) * (p.leg_len^2 + p.leg_dia^2);
p.I = (p.m/3) * (p.leg_len^2);

% Moment of water inertia (kg m^2)
p.I_water = pi/3 * p.rho * (p.leg_dia/2)^2 * p.leg_len^3;

% Length of out-lever (m)
p.hOut = sqrt(p.I/p.m);

% Drag product
p.drag_prod = (pi/4) * p.Cd * p.rho * (p.leg_dia/2) * p.leg_len^4;


%p.hOut
p.hIn/p.hOut

%% Non-dimensionalize parameters

% Constant to non-dimensionalize time
t_const = 1/sqrt(p.k/p.m);

% Specify global parameters that factor into governing equation
global I I_water hIn hOut k xJoint g drag_prod m leg_len max_theta

% Non-dimensionalize parameter values
m         = p.m/p.m;
I         = p.I/p.m/(p.hIn+p.hOut)^2;
I_water   = p.I_water/p.m/(p.hIn+p.hOut)^2;
drag_prod = p.drag_prod/p.m/(p.hIn+p.hOut)^2;
k         = p.k/p.m*t_const^2;
hIn       = p.hIn/(p.hIn+p.hOut);
hOut      = p.hOut/(p.hIn+p.hOut);
xJoint    = p.xJoint/(p.hIn+p.hOut);
g         = p.g/(p.hIn+p.hOut)*t_const^2;
leg_len   = p.leg_len/(p.hIn+p.hOut);
max_theta = p.max_theta;


%% Simulation parameters

tspan      = [0 100*10^-3]./t_const;
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

% Make sure sim ran for long enough
if t(end)==tspan(2)
    error('Evaluation time not set long enough')
end 

% Position of intersection of spring and lever
pIn  = [ p.xJoint+p.hIn*cos(theta)   p.hIn*sin(theta)];

% Position of center of intertia
pOut = [ p.xJoint-p.hOut*cos(theta)   -p.hOut*sin(theta)];

% Position of center of mass
pCOM = [ p.xJoint-p.leg_len/2*cos(theta) -p.leg_len/2*sin(theta)];

% Length of spring
%L  = pIn(:,1)-p.xSpring;
L = p.xJoint + p.hIn*cos(theta);

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

if visData
    disp(' ')
    disp(['max rotation = ' num2str(max(Dtheta).*180/pi/1000) ' deg/ms'])
    disp(' ')
    % disp(' ')
    % disp(['vmax = ' num2str(max(v))])
    % disp(['Prediction = ' num2str(sqrt(p.k/p.m)*p.L_s)])
    % disp(' ')
    iMax = find(theta==max(theta),1,'first');
    disp(['Duration of kick = ' num2str(t(iMax)*1000) ' ms'])
    disp(' ')
    disp(['Range of theta = ' num2str(range(theta)*(180/pi)) ' deg'])
    disp(' ')
end


%% Analytical approximation

%const       = sin(20/180*pi);         
% theta_ana =  p.L_s/p.hIn + (-p.L_s + p.hIn*p.theta0)/p.hIn ...
%             .* cos(sqrt(const*p.k/p.m)*(p.hIn/p.hOut).*t);
%         
% Dtheta_ana = sqrt(const*p.k/p.m) .* (p.L_s - p.hIn.*p.theta0)/p.hOut ...
%             .* sin(sqrt(const*p.k/p.m)*(p.hIn/p.hOut).*t);

hIn_const = p.hIn*sin(30/180*pi);

theta_ana =  p.L_s/hIn_const + (-p.L_s + hIn_const.*p.theta0)/hIn_const ...
            .* cos(sqrt(p.k/p.m)*(hIn_const/p.hOut).*t);
        
% Dtheta_ana = sqrt(p.k/p.m) .* (p.L_s - hIn_const.*p.theta0)/p.hOut ...
%             .* sin(sqrt(p.k/p.m)*(hIn_const/p.hOut).*t);

lambda = pi/(p.hIn/p.hOut) * sqrt(p.m/p.k);

Dtheta_ana = p.L_s.*sqrt(p.k/p.m)/p.hOut .* (0.5-0.5.*cos(pi.*t./lambda));
        

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
    plot(t,theta.*(180/pi),'-')
%     hold on
%     plot(t,theta_ana.*(180/pi),'r--')
%     hold off
    %grid on
    xlabel('Time (s)')
    ylabel('Theta (deg)')
    %legend('numerical','analytical','Location','NorthWest')
xlim([0 6e-3])
ylim([0 180])
    
    subplot(4,1,2)
    plot(t,Dtheta.*(180/pi/1000),'-')
%     hold on
%     plot(t,Dtheta_ana.*(180/pi/1000),'r--')
%     hold off
    xlabel('Time (s)')
    ylabel('dTheta/dt (deg/ms)')
    %grid on
xlim([0 6e-3])
ylim([0 90])
    
    subplot(4,1,3)
    plot(t,v,'-')
    xlabel('Time (s)')
    ylabel('Speed (m/s)')
    %grid on
xlim([0 6e-3])
ylim([0 40])
    
    subplot(4,1,4)
    plot(t,E_kinetic.*1000,'-',t,E_potential.*1000,'r',t,E_tot.*1000,'k--')
%     hold on
%     plot(t,E_g.*1000,'g',t,E_drag.*1000,'m')
%     hold off
xlim([0 6e-3])
    
    %grid on
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
d.L_range = range(L);

if d.theta(end) >= p.max_theta
    d.hit_max = true;
else
    d.hit_max = false;
end


function dX = gov_eqn(t,X)
global I I_water hIn hOut k xJoint g drag_prod m leg_len

% Angle of lever system
theta  = X(1);

% Rate of angular rotation of lever system
Dtheta = X(2);

% Position of intersection of spring and lever
pIn  = [ xJoint+hIn*cos(theta)   hIn*sin(theta)];

% Position of center of intertia
pOut = [ xJoint-hOut*cos(theta)   -hOut*sin(theta)];

% Position of center of mass
pCOM = [ xJoint-leg_len/2*cos(theta) -leg_len/2*sin(theta)];

% Current length of spring
%L  = max([pIn(1)-pSpring(1) 0]);
%L  = pIn(1)-pSpring(1);
L = xJoint + hIn*cos(theta);

% Torque generated by spring
%tau_in = hIn*Fin.*sin(theta);
tau_spring = cross([pIn 0],[-k*L 0 0]);
%tau_spring = cross([0 hIn*sin(theta) 0],[-k*L 0 0]);

% Torque from weight
tau_g = cross([pCOM 0],[0 -m*g 0]);

% Torque from drag
tau_drag = -drag_prod .* Dtheta .* abs(Dtheta);

% Define output: rate of rotation
dX(1,1) = Dtheta;

% Define output: rotational acceleration
dX(2,1) = (tau_spring(3) + tau_g(3) + tau_drag)/(I + I_water);




function [value,isterminal,direction] = evnts(t,X)

% Define global variables to be used
global xJoint hIn max_theta

% Angle of lever system
theta  = X(1);

% Rate of angular rotation of lever system
Dtheta = X(2);

% Halts execution of the model
isterminal = 1;


% Current length of spring
L = xJoint + hIn*cos(theta);

if (theta>=max_theta) || (L<=0)
    value = 0;
    direction = 0;
else
    value = 1;
    direction = 1;
end

% %if hIn*cos(X(1))<pSpring(1)
% if X(2)<=0
%     value = 0;
%     direction = 0;
% elseif X(1) >=max_theta
%     value = 0;
%     direction = 0;
% else
%     value = 1;
%     direction = 1;
% end