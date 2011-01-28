function [d,result] = run_sim(p,echoData)
% This collects parameters values, saves them to disk for the Mathematica
% kernel to find, runs the kernel, reads the resutls from disk, and passes
% them to the output structure d. 


%% Default inputs

result = {};

if nargin < 2
    echoData = 1;
end


%% Check input geometry - return empty "d", if geometry not possible

L = check_linkage(p,0);

% Report problems with the geometry
if (p.thetaStart < L.thetaInMin)
    result{length(result)+1} = ...
     'Starting theta value below what is possible for linkage geometry';

elseif (p.thetaStart > L.thetaInMax)
   result{length(result)+1} = ...
     'Starting theta value above what is possible for linkage geometry';  
end

if isempty(L)
    result{length(result)+1} = ...
     'No range of motion possible for requested geomtry';  
end

% Stop execution if there is a reported problem
if ~isempty(result)
    d = [];
    beep;beep;
    warning('No simulation run')
    disp(' ')
    disp(result)
    disp(' ')
    return
end

clear L


%% Scale input parameter values for numerical stability

% Parameter structure 's' -- scaled version of 'p'
global s

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
Lout = (sqrt((p.dacI)/p.dacMass));
s.Lout   = Lout * sL;

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
gamma_0 = atan2((B(:,1)-C(:,1)),-(B(:,2)-C(:,2)));

% Output angle
%phi_0 = acos((s.L3^2 + s.L4^2 - h_BD^2) / (2*s.L3*s.L4));

clear h_BD B C si_0


%% Run simulation

% Simulation parameters
options    = odeset('Events',@evnts,'RelTol',s.rel_tol);
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
d.theta = calc_theta(d.gamma,p.L1,p.L2,p.L3,p.L4);


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

% Angle of link 3 in gobal FOR
%eta = atan2((d.B(:,1)-d.C(:,1)),-(d.B(:,2)-d.C(:,2)));

% % Kinematic transmission
% d.KT = (p.L1*p.L2)/(p.L3*p.L4) .* csc(d.phi) .* ...
%        sqrt(1-((p.L1^2+p.L2^2-p.L3^2-p.L4^2+2*p.L3*p.L4.*cos(d.phi)).^2) ...
%        ./ (4*p.L1^2 *p.L2^2));
 
% Calculate KT for all t using a 5th-order polynomial
coef = polyfit(d.theta,d.gamma,5);
d.KT = polyval(coef,d.theta);

% Spring torque
d.tau_spring = -p.kSpring*(p.thetaRest - d.theta);

% Elastic energy
d.E_elastic  = 0.5 .* p.kSpring .* (p.thetaRest - d.theta ).^2;

% Kinetic energy
d.E_kin = 0.5 * (p.dacI) .* d.Dgamma.^2;
%d.E_kin = [0; 0.5 * (p.dacI+p.waterI) .* (diff(d.eta)./diff(d.t)).^2];

% Calculate drag energy
%d.E_drag = cumtrapz(abs(d.phi(1)-d.phi),abs(d.dragTau));
d.E_drag = zeros(size(d.t,1),size(d.t,2));   



          
          
function dX = gov_eqn(t,X)
global s

% Output angle of lever system
gamma  = X(1);

% Rate of rotation of input angle
Dgamma = X(2);

% Calculate theta for current gamma
theta = calc_theta(gamma,s.L1,s.L2,s.L3,s.L4);

% Length between points b & d 
h_BD = sqrt(s.L1^2 + s.L2^2 - 2*s.L1*s.L2.*cos(theta));

% Angle btwn links 1 & 4
si = acos((h_BD.^2 + s.L1^2 - s.L2^2)./(2*h_BD.*s.L1)) + ...
     acos((h_BD.^2 + s.L4^2 - s.L3^2)./(2*h_BD.*s.L4));
   
% Position vector for point B
B = [s.L2.*sin(theta) s.L2.*cos(theta)];

% Position vector for point C
C = [s.L4.*sin(si) s.L1-s.L4.*cos(si)];

% Torque created by spring
tau_spring = -s.kSpring*(s.thetaRest - theta);

% Unit vector for input force at B (perpendicular to link 2)
F_unit(1,1) =  B(2)./ sqrt(B(1)^2 + B(2)^2);
F_unit(1,2) = -B(1)./ sqrt(B(1)^2 + B(2)^2);

% Spring force vector created at position B
F_B(1,1) = (tau_spring/s.L2) * F_unit(1);
F_B(1,2) = (tau_spring/s.L2) * F_unit(2);
F_B(1,3) = 0;

% Position vector for in-lever arm
L_B(1,1) = C(1) - B(1);
L_B(1,2) = C(2) - B(2);
L_B(1,3) = 0;

dacI = (s.dMass .* (s.Lout)^2);

% Torque driving the motion of the carpus
tau_in = cross(L_B,F_B);

% Collapse dimensions of the torque vector
tau_in = tau_in(3);

% Define output: rate of rotation
dX(1,1) = Dgamma;

% Define output: rotational acceleration
dX(2,1) = tau_in/(s.dacI);



function [value,isterminal,direction] = evnts(t,X)

% Define global variables to be used
global s

% Output angle of lever system
gamma  = X(1);

% Rate of rotation of input angle
Dgamma = X(2);

% Calculate theta for current gamma
theta = calc_theta(gamma,s.L1,s.L2,s.L3,s.L4);

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
    error('Sim stopped early: h_BD cannot exceed L3 + L4')
elseif h_BD < abs(s.L4-s.L3)
    value = 0;
    direction = 0;
    error('Sim stopped early: h_BD cannot be less than L4 - L3')       
elseif h_BD > (s.L1 + s.L2)
    value = 0;
    direction = 0;
    error('Sim stopped early: h_BD cannot exceed L1 + L2')
elseif h_BD < abs(s.L1 - s.L2)
    value = 0;
    direction = 0;
    error('Sim stopped early: h_BD cannot be less than L1 - L2')
elseif (theta>=s.thetaRest)
    value = 0;
    direction = 0;
else
    value = 1;
    direction = 1;
end


function theta = calc_theta(gamma,L1,L2,L3,L4)
theta = acos((L1.^3.*L2 + L1.*L2.^3 + 2.*L1.*L2.*L3.^2 - L1.*L2.*L4.^2 - ... 
        L2.*L3.*(3.*L1.^2 + L2.^2 + L3.^2 - L4.^2).*cos(gamma) + ...
        L1.*L2.*L3.^2.*cos(2.*gamma) + ...
        sqrt(-L2.^2.*L3.^2.*(L1.^4 - 2.*L1.^2.*L2.^2 + L2.^4 + ...
         4.*L1.^2.*L3.^2 - 2.*L2.^2.*L3.^2 + L3.^4 - 2.*L1.^2.*L4.^2 - ... 
         2.*L2.^2.*L4.^2 - 2.*L3.^2.*L4.^2 + L4.^4 - ...
         4.*L1.*L3.*(L1.^2 - L2.^2 + L3.^2 - L4.^2).*cos(gamma) + ... 
         2.*L1.^2.*L3.^2.*cos(2.*gamma)).*sin(gamma).^2))./ ...
         (2.*L2.^2.*(L1.^2 + L3.^2 - 2.*L1.*L3.*cos(gamma))));
     %theta = pi/2-theta;
