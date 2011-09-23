function d = run_simple
% Runs a numerical simulation of the predatory strike of a mantis shrimp. 
%
% This code was developed for the McHenry, Claverie, Rosario & Patek (in
% review) paper.  Parameter values correspond to Individual #2 reported in
% that study.


%% Simulation parameters

% Whether to display details of a simulation
echoData = 1;

% Relative tolerence of the simulation
p.rel_tol = 1e-7;

% Precision of simulation (default = 10^-5, use 10^-7 for high precision)
p.maxError = 10^-7;

% Maximum duration of simulation (s)
p.simDur    = 0.02;

% Maximum step size of simulation (s)
p.maxStep   = 1e-5;


%% Morphological and mechanical parameters

% Density of water (kg/m^3)
p.rho = 998;

% Below are values specific to Individual #2

% Linkage lengths (m)
p.L1 = 6.536e-3;
p.L2 = 4.227e-3;
p.L3 = 0.5676e-3;
p.L4 = 7.5588e-3;
        
% Dactyl length (m) -- calculated as the length btwn point b and the
% distal end of the dactyl
p.dac_len = 10.902e-3;   

% Mass of the striking body (kg) -- estimated from scaling relationships of
% CT specimens
p.dac_mass =  10^-6 .* 10.^(3.13.*log10(p.dac_len.*1000) -0.885);           
           
% Moment of inertia for striking body (kg m^2) 
% (calculated using I* = 0.28608)
p.dac_I = 0.28608 .* ((p.dac_len).^2) .* p.dac_mass;

% Initial output angle of the dactyl (rad)
% Equal to the angle btwn link 4 and axis btwn f and g
p.gamma = (pi/180).*60.0725;

% Initial input angle (rad)
p.thetaStart = (pi/180).*80;

% Resting angle for torsion spring (rad)
p.thetaRest = (pi/180).*92.31;
                   
% Linear spring stiffness (N/m)               
p.k_lin = 27.322e3;
             
% Angle between link 2 and wire during materials testing (rad)
lambda = (pi/180).*25.787;  

% Distance between load and point A (m)
dist_load = 4.821e-3;

% Torsion spring at mV joint (Nm/rad) 
p.k_spring = p.k_lin .* dist_load.^2 .* cos(lambda);

% Dimensionless drag index
p.D = 0.0756;

% Moment of inertia for the water
p.waterI = 0.0539e-8;         


% Clear variables not used below
clear lambda dist_load


%% Global variables declared
% These variables are passed to the governing function during the
% simulation

global s sp_theta sp_KT


%% Check input geometry & calculate splines

[sp_gamma,sp_KT,sp_theta,min_theta,max_theta] = ...
                                    linkage_functions(p.L1,p.L2,p.L3,p.L4);

if p.thetaStart < min_theta
    error('Starting value for theta not theoretically possible');
    
elseif p.thetaRest > max_theta
    error('Resting theta value greater than what is theoretically possible');

end


%% Scale input parameter values for numerical stability
% All parameters used by the model are rescaled, made dimensionless, and
% stored in the 's' structure.

% Scaling factors 
sL = 1 / p.L1;
sM = 1/p.dac_mass;
sT = 10^3;

% Dimensionless parameters
s.gamma       = p.gamma;
s.thetaStart  = p.thetaStart;
s.thetaRest   = p.thetaRest;
s.DrgIdx      = p.D;
s.rel_tol     = p.rel_tol;

% Linear dimensions
s.L1     = p.L1 * sL;
s.L2     = p.L2 * sL;
s.L3     = p.L3 * sL;
s.L4     = p.L4 * sL;

% Mechanical properties
s.dMass   = p.dac_mass * sM;
s.dacI    = p.dac_I * sM * sL^2;
s.waterI  = p.waterI * sM * sL^2;
s.kSpring = p.k_spring * (sM * sL^2/sT^2);
s.rho     = p.rho * sM / sL^3;

% Time
s.simDuration = p.simDur * sT;
s.maxStep     = p.maxStep * sT;


%% Define initial geometry
% Necessary for the initial conditions of the simulation

% Calculate length between B & D for range of theta
h_BD     = sqrt(s.L1^2 + s.L2^2 - 2.*s.L1.*s.L2.*cos(s.thetaStart));

% Check that the initial linkage geometry is possible
if h_BD > (s.L3 + s.L4)
    error('h_BD cannot exceed L3 + L4')
elseif h_BD < abs(s.L4-s.L3)
    error('h_BD cannot be less than L4 - L3')       
elseif h_BD > (s.L1 + s.L2)
    error('h_BD cannot exceed L1 + L2')
elseif h_BD < abs(s.L1 - s.L2)
    error('h_BD cannot be less than L1 - L2')
end

% Angle btwn links 1 & 4
si_0 = acos((h_BD.^2 + s.L1^2 - s.L2^2)./(2*h_BD.*s.L1)) + ...
       acos((h_BD.^2 + s.L4^2 - s.L3^2)./(2*h_BD.*s.L4));

% Positions of points B & C
B = [s.L2.*sin(s.thetaStart) s.L2.*cos(s.thetaStart)];
C = [s.L4.*sin(si_0) s.L1-s.L4.*cos(si_0)];

% Initial value for output angle
gamma_0 = 3*pi/2 + atan2((B(:,1)-C(:,1)),-(B(:,2)-C(:,2)));

% Clear extra variables
clear h_BD B C si_0


%% Run simulation

% Simulation parameters
options    = odeset('Events',@evnts,'RelTol',s.rel_tol,...
                    'MaxStep',s.maxStep);
tspan      = [0 s.simDuration];
init_vals  = [gamma_0; 0];

% Intialize timer
if echoData
    tic;
end

% Run
[t,X] = ode45(@gov_eqn,tspan,init_vals,options);

% Update sim time
if echoData
    tlapse = toc;
    disp(['Sim complete in ' num2str(tlapse) ' s'])
end

% Store results
d.t      = t ./ sT;
d.gamma  = X(:,1);
d.Dgamma = X(:,2) .* sT;

% Clear others
clear t X tspan init_vals s sT sL sM


%% Calculate result variables

% Calculate theta for gamma values
d.theta = fnval(sp_theta,d.gamma);

% Length between points b & d 
h_BD = sqrt(p.L1^2 + p.L2^2 - 2*p.L1*p.L2.*cos(d.theta));

% Angle btwn links 1 & 4
d.si = acos((h_BD.^2 + p.L1^2 - p.L2^2)./(2*h_BD.*p.L1)) + ...
       acos((h_BD.^2 + p.L4^2 - p.L3^2)./(2*h_BD.*p.L4));

% Positions of points
d.A = zeros(length(d.t),2);
d.B = [p.L2.*sin(d.theta) p.L2.*cos(d.theta)];
d.C = [p.L4.*sin(d.si) p.L1-p.L4.*cos(d.si)];
d.D = [zeros(length(d.t),1) p.L1.*ones(length(d.t),1)];

% Find KT for theta values
d.KT = fnval(sp_KT,d.theta);

% Spring torque
d.tau_spring = p.k_spring*(p.thetaRest - d.theta);

% Drag torque
d.tau_drag = 0.5*p.rho.*d.Dgamma.^2 .* p.dac_len^5 .* p.D  .* ... 
             d.Dgamma./abs(10^-20 + d.Dgamma);

% Elastic energy
d.E_elastic  = 0.5 .* p.k_spring .* (p.thetaRest - d.theta ).^2;

% Kinetic energy
d.E_kin = 0.5 * (p.dac_I+p.waterI) .* d.Dgamma.^2;

% Drag energy
d.E_drag = cumtrapz(abs(d.gamma-d.gamma(1)),abs(d.tau_drag));

% Striking body momentum
d.p_SB = d.Dgamma.*p.dac_I/sqrt(p.dac_I/p.dac_mass);
          

%% Plot simulation results
figure;

% Input and output angles
subplot(4,1,1)
plot(d.t.*1000,(d.gamma)*180/pi,'b',...
     d.t.*1000,d.theta.*180/pi,'r');
ylabel('angle (deg)')
xlabel('time (ms)')

legend('gamma','theta')
grid on

% Energy
subplot(4,1,2)
E_tot = d.E_kin + d.E_drag + d.E_elastic;
plot(d.t.*1000,d.E_elastic.*1000,'r',...
     d.t.*1000,d.E_drag.*1000,'b',...
     d.t.*1000,d.E_kin.*1000,'g',...
     d.t.*1000,E_tot.*1000,'k--')
ylabel('Energy (mJ)')
xlabel('time (ms)')
grid on
legend('elastic','drag','kinetic','total','Location','West');
title(['min KT = ' num2str(min(d.KT))])

% Power 
subplot(4,1,3)
plot(d.t(2:end).*1000,diff(d.E_elastic)./diff(d.t),'r',...
    d.t(2:end).*1000,diff(d.E_drag)./diff(d.t),'b',...
    d.t(2:end).*1000,diff(d.E_kin)./diff(d.t),'g')
ylabel('Power (J/s)')
xlabel('time (ms)')
grid on

% Momentum
subplot(4,1,4)
plot(d.t.*1000,d.p_SB.*1000,'k')
grid on
xlabel('time (ms)')
ylabel('Striking body momentum (g m/s)')




%% FUNCTIONS
          
function dX = gov_eqn(t,X)
% Defines the differential equations used to calculate the instananeous
% torques acting on the striking body during a strike

% Global variables allow parameter values to be passed into this function
global s sp_theta sp_KT

% Output angle of lever system
gamma  = X(1);

% Rate of rotation of output angle
Dgamma = X(2);

% Use previously calculated spline to determine theta for current gamma
theta = fnval(sp_theta,gamma);

% KT calculated from spline for current theta
KT = fnval(sp_KT,theta);

% Torque created by spring
tau_spring = -s.kSpring*(s.thetaRest - theta);

% Torque created by drag
tau_drag = -0.5*s.rho.*Dgamma.^2 .* (s.dacLen^5) .* s.DrgIdx  .* ... 
             Dgamma./abs(10^-20 + Dgamma);
         
% Input torque
tau_in = - (tau_spring./KT);

% Define output: rate of rotation
dX(1,1) = Dgamma;

% Define output: rotational acceleration
dX(2,1) = (tau_in + tau_drag) ./ (s.dacI+s.waterI);


function [value,isterminal,direction] = evnts(t,X)
% Checks the status of a simulation to ensures that the geometry is
% possible and halts a simulation when the resting length of the spring is
% reached.

% Define global variables to be used
global s sp_theta

% Output angle of lever system
gamma  = X(1);

% Rate of rotation of input angle
Dgamma = X(2);

% Calculate theta for current gamma
theta = fnval(sp_theta,gamma);
%theta = calc_theta(gamma,s.L1,s.L2,s.L3,s.L4);

% Halts execution of the model
isterminal = 1;

% Length between points b & d 
h_BD = sqrt(s.L1^2 + s.L2^2 - 2*s.L1*s.L2.*cos(theta));

% Input angle
theta = acos((s.L1^2 + s.L2^2 - h_BD^2) / (2*s.L1*s.L2));

% Check that linkage geometry is possible
if h_BD > (s.L3 + s.L4)
    value = 0;
    direction = 0;
    error('Sim stopped: h_BD cannot exceed L3 + L4')
elseif h_BD < abs(s.L4-s.L3)
    value = 0;
    direction = 0;
    error('Sim stopped: h_BD cannot be less than L4 - L3')       
elseif h_BD > (s.L1 + s.L2)
    value = 0;
    direction = 0;
    error('Sim stopped: h_BD cannot exceed L1 + L2')
elseif h_BD < abs(s.L1 - s.L2)
    value = 0;
    direction = 0;
    error('Sim stopped: h_BD cannot be less than L1 - L2')
elseif (theta>=s.thetaRest)
    value = 0;
    direction = 0;
else
    value = 1;
    direction = 1;
end


function [sp_gamma,sp_KT,sp_theta,min_theta,max_theta] = linkage_functions(L1,L2,L3,L4)
% Returns spline functions that describe how gamma and KT vary with theta.
% Also calculates the range of theoretically possible theta values.

% Number of discrete values to consider
numVals = 10000;

% Define full range of possible values for gamma
theta = linspace(0,pi,numVals)';

% Calculate length between B & D for range of theta
h_BD     = sqrt(L1^2 + L2^2 - 2.*L1.*L2.*cos(theta));

% Check range of possible gamma
idx  = (h_BD < (L3 + L4)) & ...
       (h_BD > abs(L4-L3)) & ...
       (h_BD < (L1 + L2)) & ...
       (h_BD > abs(L1 - L2));

% Quit, if no range of motion is possible
if max(idx) == 0  
    error('No input angles possible with this geometry');
end

% Trim impossible range of motion values
theta = theta(idx);
h_BD  = h_BD(idx);

clear idx

% phi - the angle between links 3 and 4
phi  = acos( (L3.^2 + L4.^2 - h_BD.^2) ./ (2.*L3.*L4) );

% psi - the angle between links 1 and 4
psi = acos((h_BD.^2 + L1^2 - L2^2)./(2*h_BD.*L1)) + ...
       acos((h_BD.^2 + L4^2 - L3^2)./(2*h_BD.*L4));

% Positions of points
A = zeros(length(theta),2);
B = [L2.*sin(theta) L2.*cos(theta)];
C = [L4.*sin(psi) L1-L4.*cos(psi)];
D = [zeros(length(theta),1) L1.*ones(length(theta),1)];

% Numerical calculation of gamma
gamma = 3*pi/2 + atan2((B(:,1)-C(:,1)),-(B(:,2)-C(:,2)));

% Spline fit to gamma values
sp_gamma = csapi(theta,gamma);

% Spline fit to theta values
sp_theta = csapi(gamma,theta);

% Find KT, the first derivative of gamma spline
sp_KT = fnder(sp_gamma);

% Range of possible theta values
min_theta = min(theta);
max_theta = max(theta);

             