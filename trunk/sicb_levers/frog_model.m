function r = frog_model(MA,b_eff)
% Model based on Lappin et al. (2006) 


%% Parameter values for the frog

% Effective mass (kg)
m_eff = 10.33/1000;

% Spring constant (N/m)
k_t = 581;

% Initial muscle length (m)
x0 = 3.2e-3./100; % High effort
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
    plot_results = 1;
    
end


%% Simulation parameters

% Relative tolerence of the simulation
rel_tol = 1e-7;

% Precision of simulation (default = 10^-5, use 10^-7 for high precision)
maxError = 10^-7;

% Maximum step size of simulation (s)
maxStep   = 1e-5;

% Maximum duration of simulation (s)
t_max     = 40e-3;


%% Calculated parameters

% Calculate mass of jaw
% The paper included the mass of the muscle in its calculation of effective
% mass, but did not report individual values for the jaw and muscle.  I
% have therefore assumed that the jaw mass is much larger than the
% effective mass.
m_jaw = m_eff * (1/8);


%% Run simulation

% Simulation parameters
options    = odeset('Events',@evnts,'RelTol',rel_tol,'MaxStep',maxStep);
tspan      = [0 t_max];
init_vals  = [x0; 0];

% Define global parameter values
global p

p.m_eff = m_jaw./MA;
p.k_t   = k_t;
p.MA    = MA;
p.b_eff = b_eff;

clear MA k_t m_jaw b_eff m_eff 

% Run
[t,X] = ode45(@gov_eqn,tspan,init_vals,options);

% Store time
r.t = t;

% Store total muscle position and velocity
r.x_t = X(:,1);
r.v_t = X(:,2);

% Calculate acceleration
r.a_t = -(p.b_eff.*r.v_t + p.k_t.*r.x_t)./p.m_eff;

% Jaw velocity
r.v_jaw = r.v_t ./ p.MA;

% Jaw acceleration
r.a_jaw = r.a_t ./ p.MA;

clear X t

% Dashpot force
r.F_dash = r.v_t.*p.b_eff;

% Spring force
r.F_spring = r.x_t.*p.k_t;

% Elastic energy
r.E_elastic  = 0.5 .* p.k_t .* r.x_t.^2;

% Kinetic energy
r.E_kin = 0.5*p.m_eff.*r.v_t.^2;

% Dashpot energy
r.E_dash = cumtrapz(abs(r.x_t),r.F_dash);

% Total energy
r.E_tot = r.E_elastic + r.E_kin + r.E_dash;

% Efficiency
r.eff = max(r.E_kin)./max(r.E_elastic);

% Analytical estimate of max velocity
r.v_t_max = sqrt(2*max(r.E_elastic)/p.m_eff);


%% Plot results


if plot_results
    
    figure
    subplot(4,1,1)
    plot(r.t.*1000,r.x_t.*1000,'k')
    xlabel('time (ms)')
    ylabel('Muscle position (mm)')
    
    subplot(4,1,2)
    plot(r.t.*1000,r.v_t.*1000,'k')
    %hold on;plot(r.t.*1000,r.v_jaw.*1000,'r')
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
global p

% Muscle position
x  = X(1);

% Contraction speed
v = X(2);

% Define output: contraction speed
dX(1,1) = v;

% Define output: contraction acceleration
dX(2,1) = -(p.b_eff.*v + p.k_t.*x)./p.m_eff;



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



