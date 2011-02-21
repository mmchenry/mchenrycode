function [d,result] = run_sim(p,echoData)
% This collects parameters values, saves them to disk for the Mathematica
% kernel to find, runs the kernel, reads the resutls from disk, and passes
% them to the output structure d. 


%% Default inputs

result = {};

if nargin < 2
    echoData = 1;
end

%% Global variables declared
% Parameter structure 's' -- scaled version of 'p'
global s sp_theta sp_KT


%% Check input geometry - return empty "d", if geometry not possible

[sp_gamma,sp_KT,sp_theta,min_theta,max_theta] = ...
                                    linkage_functions(p.L1,p.L2,p.L3,p.L4);

if p.thetaStart < min_theta
    error('Starting value for theta not theoretically possible');
    
elseif p.thetaRest > max_theta
    error('Resting theta value greater than what is theoretically possible');

end


%% Scale input parameter values for numerical stability

% Scaling factors 
sL = 1 / p.L1;
sM = 1/p.dacMass;
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
s.L5     = p.L5 * sL;
s.hAF    = p.h_AF * sL;
s.dacLen = p.dacLen * sL;
%Lout = (sqrt((p.dacI)/p.dacMass));
%s.Lout   = Lout * sL;

% Mechanical properties
s.dMass   = p.dacMass * sM;
s.dacI    = p.dacI * sM * sL^2;
%s.dacI    = (p.dacMass .* (Lout)^2) * sM * sL^2;
s.waterI  = p.waterI * sM * sL^2;
s.kSpring = p.kSpring * (sM * sL^2/sT^2);
s.rho     = p.rho * sM / sL^3;

% Time
s.simDuration = p.simDur * sT;


%% Define initial geometry

% Calculate length between B & D for range of theta
h_BD     = sqrt(s.L1^2 + s.L2^2 - 2.*s.L1.*s.L2.*cos(s.thetaStart));

% Check that linkage geometry is possible
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

% Initial value for eta
%gamma_0 = atan2((B(:,1)-C(:,1)),-(B(:,2)-C(:,2)));
gamma_0 = 3*pi/2 + atan2((B(:,1)-C(:,1)),-(B(:,2)-C(:,2)));

% Output angle
%phi_0 = acos((s.L3^2 + s.L4^2 - h_BD^2) / (2*s.L3*s.L4));

clear h_BD B C si_0


%% Run simulation

% Simulation parameters
options    = odeset('Events',@evnts,'RelTol',s.rel_tol,...
                    'MaxStep',s.simDuration./1000);
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


%% Calculate result parameters

% Calculate theta for gamma values
%d.theta = calc_theta(d.gamma,p.L1,p.L2,p.L3,p.L4);
d.theta = fnval(sp_theta,d.gamma);

% % Rate of change in input angle
% d.Dtheta = d.Dphi .* (p.L3 * p.L4 .* sin(d.phi)) ./ ...
%          (p.L1.*p.L2.*sqrt(1 - ...
%          ((p.L1^2 + p.L2^2 - p.L3^2 - p.L4^2 + 2*p.L3*p.L4 .*cos(d.phi)).^2)./...
%          (4*p.L1^2*p.L2^2)));

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

% Find KT for theta value
KT = fnval(sp_KT,d.theta);

% Calculate KT for all t using a 5th-order polynomial
cs = fnder(csapi(d.theta,d.gamma));
d.KT = fnval(cs,d.theta);
clear cs

% Spring torque
d.tau_spring = p.kSpring*(p.thetaRest - d.theta);

% Drag torque
d.tau_drag = 0.5*p.rho.*d.Dgamma.^2 .* p.dacLen^5 .* p.D  .* ... 
             d.Dgamma./abs(10^-20 + d.Dgamma);

% Elastic energy
d.E_elastic  = 0.5 .* p.kSpring .* (p.thetaRest - d.theta ).^2;

% Kinetic energy
d.E_kin = 0.5 * (p.dacI+p.waterI) .* d.Dgamma.^2;
%d.E_kin = [0; 0.5 * (p.dacI+p.waterI) .* (diff(d.eta)./diff(d.t)).^2];

% Calculate drag energy
d.E_drag = cumtrapz(abs(d.gamma-d.gamma(1)),abs(d.tau_drag));
%d.E_drag = zeros(size(d.t,1),size(d.t,2));   
          
          
function dX = gov_eqn(t,X)
global s sp_theta sp_KT

% Output angle of lever system
gamma  = X(1);

% Rate of rotation of input angle
Dgamma = X(2);

% Calculate theta for current gamma
%theta = calc_theta(gamma,s.L1,s.L2,s.L3,s.L4);
theta = fnval(sp_theta,gamma);

% KT
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




             