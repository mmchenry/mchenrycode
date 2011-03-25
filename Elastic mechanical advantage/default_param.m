function p = default_param
% Default parameter values used in the locust model

% Parmeter values ('leg' = tibia + tarsus)
%   leg mass = 24e-6 (Bennet-Clark, 1975)
%   leg diameter = 1.18e-3; (Bennet-Clark, 1975)
%   leg length = 30e-3; (Bennet-Clark, 1975)
%   k spring (estimate) = 1.2e4; (Bennet-Clark, 1975)
%   Shortening length = 1.02e-3; (Bennet-Clark, 1975)
%   In-lever length = 0.76e-3 (Heitler, 1974)
%   Initial angle btwn tendon & leg = 6 deg (Heitler, 1974)

% Mass (kg)
p.m = 21e-6;

% Spring stiffness (N/m)
p.k = 0.3e4;

% Spring stretch -- initial spring length beyond resting length (m)
p.L_s = 1.5e-3; % Makes and averge kick with  p.k = 0.3e4;
%p.L_s = 0.97e-3; % Makes an 'average' kick with k = 0.73e4;
%p.L_s = 1.2e-3; % Makes a 3 ms kick with k = 2e4

% Distance btwn spring insertion and pivot (m)
p.hIn = 0.76e-3;

% Diameter of leg (m)
p.leg_dia = 1.18e-3;

% Length of leg (m)
p.leg_len = 30e-3;

% Initial angle of lever system (rad)
p.theta0 = 6/180*pi;

% Acceleration of gravity (m/s^2)
p.g = 0*9.81;

% Density of air (kg m^-3)
p.rho = 0*1.2;

% Dynamic viscosity of air (Pa s)
p.mu = 1.845e-5;

% Drag coefficient for a cylinder (dimensionless)
p.Cd = 1.2;

% Max excursion of theta (rad)
p.max_theta = 200/180*pi;
