function frog_model_play
% Model based on Lappin et al. (2006) 


% Current individual
%indiv = 1;
%[m_eff,MA,k,b_eff] = get_param(indiv);

% Initial position
%x0 = 3.5e-3;

% Global variables allow parameter values to be passed into this function
global b_eff k_t m_eff


% Relative tolerence of the simulation
s.rel_tol = 1e-7;

% Precision of simulation (default = 10^-5, use 10^-7 for high precision)
s.maxError = 10^-7;

% Maximum step size of simulation (s)
s.maxStep   = 1e-5;


[m_eff,MA,k_t,b_eff,x0] = get_param_mean;


% Calculate mass of jaw
% The paper included the mass of the muscle in its calculation of effective
% mass, but did not report individual values for the jaw and muscle.  I
% have therefore assumed that the jaw mass is much larger than the
% effective mass.
m_jaw = m_eff * MA;

% Parameters necessary to solve the kinematics
omega_n   = sqrt(k_t/m_eff);
zeta      = b_eff/(2*m_eff*omega_n);
omega_d   = sqrt(1-zeta^2*omega_n);

% Defien time vector
t = linspace(0,40e-3,10000);

% Zero out damping (for energy conservation sims)
%b_eff = 0;

% Changes in position of the muscle
x_t = exp(-zeta*omega_n.*t) .* ...
         (x0.*cos(omega_d.*t) + (omega_n*x0/omega_d).* sin(omega_d.*t));

% Change in muscle velocity over time
u_t = (-exp(-zeta*omega_n.*t)./omega_d) .* ...
         (x0.*((zeta-1).*omega_d.*omega_n.*cos(t.*omega_d) + ...
          (omega_d.^2+zeta.*omega_n.^2).*sin(t.*omega_d)));

% Change in muscle acceleration over time
a_t = (exp(-zeta*omega_n.*t)./omega_d) .* ...
          (x0.*(-omega_d.*(omega_d.^2-(zeta-2).*zeta.*omega_n.^2).*cos(t.*omega_d) + ...
                omega_n.*((2.*zeta-1).*omega_d.^2+zeta.^2.*omega_n.^2).*sin(t.*omega_d)));
      
% Jaw velocity
u_jaw = u_t ./ MA;

% Jaw acceleration
a_jaw = a_t ./ MA;



      

%% Run simulation

% Simulation parameters
options    = odeset('RelTol',s.rel_tol,...
                    'MaxStep',s.maxStep);
tspan      = [0 max(t)];
init_vals  = [x0; 0];


% Run
[t_sim,X] = ode45(@gov_eqn,tspan,init_vals,options);

x_sim = X(:,1);
v_sim = X(:,2);


% Dashpot force
F_dash = v_sim.*b_eff;

% Elastic energy
E_elastic  = 0.5 .* k_t .* x_sim.^2;

% Kinetic energy
E_kin = 0.5*m_eff.*v_sim.^2;

% Dashpot energy
E_dash = cumtrapz(abs(x_sim),F_dash);

% Total energy
E_tot = E_elastic + E_kin + E_dash;


figure

subplot(4,1,1)
plot(t.*1000,x_t.*1000,'k')
hold on; plot(t_sim.*1000,x_sim.*1000,'r--')
xlabel('time (ms)')
ylabel('Muscle position (mm)')
title('Kinematics')
legend('analytical','numerical')

u_diff = diff(x_t)./diff(t);

subplot(4,1,2)
plot(t.*1000,u_t.*1000,'k')
hold on; plot(t_sim.*1000,v_sim.*1000,'r--')
%hold on; plot(t(2:end).*1000,u_diff.*1000,'r--')
xlabel('time (ms)')
ylabel('Muscle vel (mm/s)')

a_diff = diff(u_diff)./mean(diff(t));

subplot(4,1,3)
plot(t.*1000,a_t.*1000,'k')
%hold on; plot(t(2:end-1).*1000,a_diff.*1000,'r--')
hold on;plot(t_sim(2:end).*1000,diff(v_sim)./diff(t_sim).*1000,'r--')
xlabel('time (ms)')
ylabel('Muscle accel (mm/s)')

subplot(4,1,4)
plot(t_sim.*1000,E_elastic.*1000,'k',t_sim.*1000,E_kin.*1000,'b',...
     t_sim.*1000,E_dash.*1000,'r',t_sim.*1000,E_tot.*1000,'k--')
xlabel('time (ms)')
ylabel('Energy (mJ)')
legend('elastic','kinetic','dashpot','total')
title('Energetics of numerical model')




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


function [m_eff,MA,k_t,b_eff,x0] = get_param_mean

% Effective mass (kg)
m_eff = 10.33/1000;

% Mechanical advantage
MA = 1/8;

% Spring constant (N/m)
k_t = 581;

% Damping coefficient (Ns/m)
b_eff = 4.6;

% Initial position (m)
x0 = 3.2e-3; % High effort
%x0 = 2.3e-3; % Low effort


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






