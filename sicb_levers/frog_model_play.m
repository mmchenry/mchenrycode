function r = frog_model_play(MA,b)
% Model based on Lappin et al. (2006) 


%% Parameter values

% Global variables allow parameter values to be passed into this function
global b_eff k_t m_eff

% Relative tolerence of the simulation
rel_tol = 1e-7;

% Precision of simulation (default = 10^-5, use 10^-7 for high precision)
maxError = 10^-7;

% Maximum step size of simulation (s)
maxStep   = 1e-5;

% Maximum duration of simulation (s)
t_max     = 40e-3;

% Effective mass (kg)
m_eff = 10.33/1000;

% Spring constant (N/m)
k_t = 581;

% Initial position (m)
x0 = 3.2e-3; % High effort
%x0 = 2.3e-3; % Low effort

% Assign defaults, if not provided
if nargin < 2
    
    % Damping coefficient (Ns/m)
    b_eff = 4.6;
    
    % Plot results
    plot_results = 1;
    
    if nargin < 1
        % Mechanical advantage
        MA = 1/8;
    end
else
    
    % Plot results
    plot_results = 0;
    
    % Define damping coefficient
    b_eff = b;
end


%% Calculated parameters

% Calculate mass of jaw
% The paper included the mass of the muscle in its calculation of effective
% mass, but did not report individual values for the jaw and muscle.  I
% have therefore assumed that the jaw mass is much larger than the
% effective mass.
m_jaw = m_eff * (1/8)^2;
%TODO:Revisit this calculation


%% Run simulation

% Simulation parameters
options    = odeset('Events',@evnts,'RelTol',rel_tol,'MaxStep',maxStep);
tspan      = [0 t_max];
init_vals  = [x0; 0];

% Run
[t,X] = ode45(@gov_eqn,tspan,init_vals,options);

% Store time
r.t = t;

% Extract total muscle position and velocity
r.x_t = X(:,1);
r.v_t = X(:,2);

% Solve for acceleration
r.a_t = -(b_eff.*r.v_t + k_t.*r.x_t)./m_eff;

% Jaw velocity
r.v_jaw = r.v_t ./ MA;

% Jaw acceleration
r.a_jaw = r.a_t ./ MA;

clear X t

% Dashpot force
r.F_dash = r.v_t.*b_eff;

% Spring force
r.F_spring = r.x_t.*k_t;

% Elastic energy
r.E_elastic  = 0.5 .* k_t .* r.x_t.^2;

% Kinetic energy
r.E_kin = 0.5*m_jaw.*r.v_jaw.^2;

% Dashpot energy
r.E_dash = cumtrapz(abs(r.x_t),r.F_dash);

% Total energy
r.E_tot = r.E_elastic + r.E_kin + r.E_dash;

% Efficiency
r.eff = max(r.E_kin)./max(r.E_elastic);


%% Plot results


if plot_results

subplot(4,1,1)
plot(r.t.*1000,r.x_t.*1000,'k')
xlabel('time (ms)')
ylabel('Muscle position (mm)')

subplot(4,1,2)
plot(r.t.*1000,r.v_t.*1000,'k')
hold on;plot(r.t.*1000,r.v_jaw.*1000,'r')
xlabel('time (ms)')
ylabel('Muscle vel (mm/s)')

subplot(4,1,3)
plot(r.t.*1000,r.a_t.*1000,'k')
xlabel('time (ms)')
ylabel('Muscle accel (mm/s)')

subplot(4,1,4)
plot(r.t.*1000,r.E_elastic.*1000,'k',r.t.*1000,r.E_kin.*1000,'b',...
     r.t.*1000,r.E_dash.*1000,'r',r.t.*1000,r.E_tot.*1000,'k--')
xlabel('time (ms)')
ylabel('Energy (mJ)')
legend('elastic','kinetic','dashpot','total')
title(['Efficiency = ' num2str(round(r.eff.*100)) ' %'])

end


function dX = gov_eqn(t,X)
% Defines the differential equations used to calculate the instananeous
% torques acting on the striking body during a strike

% Global variables allow parameter values to be passed into this function
global b_eff k_t m_eff

% Muscle position
x  = X(1);

% Contraction speed
v = X(2);

% Define output: contraction speed
dX(1,1) = v;

% Define output: contraction acceleration
dX(2,1) = -(b_eff.*v + k_t.*x)./m_eff;



function [value,isterminal,direction] = evnts(t,X)
% Checks the status of a simulation to ensures that the geometry is
% possible and halts a simulation when the resting length of the spring is
% reached.

% Define global variables to be used
%global b_eff k_t m_eff

% Muscle position
x  = X(1);

% Contraction speed
v = X(2);

% Halts execution of the model
isterminal = 1;

% Halt when spring at zero
if x < .001e-3
    value = 0;
    direction = 0;
else
    value = 1;
    direction = 1;
end


function [m_eff,MA,k_t,b_eff] = get_param(indiv)

% Effective mass (kg)
m_eff_tmp = [9.12 9.94 9.33 12.9]./1000;

% Mechanical advantage
MA_tmp = [1/7.9 1/7.9 1/8.6 1/7.7];

% Spring constant (N/m)
k_tmp = [712 582 522 508];

% Damping coefficient (Ns/m)
b_tmp = [5.6 4.4 4.1 4.4];

% Keep value for requested individual
m_eff = m_eff_tmp(indiv);
MA    = MA_tmp(indiv);
k_t   = k_tmp(indiv);
b_eff     = b_tmp(indiv);
