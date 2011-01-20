function run_single_sim
% Runs the model and visualizes the resutls for a single individual
% Similar to model_test, but doesn't rerun for all indivduals

% Indivudal to simulate
indiv = 120;


%% Path definitions (computer specific)

% Path to text file of Mathematica commands
p.mathFile = ['"<</Volumes/data_commuter/Projects/Patek_project/m_files_podmodel/sim_code.txt"'];

% Path to the Mathematica kernel
p.kernelPath = '/Applications/Mathematica.app/Contents/MacOS/MathKernel';

% Root path for simulation data
p.simsPath = '/Volumes/data_commuter/Projects/Patek_project/sims';

% Path to force data
forcePath = '/Volumes/data_commuter/Projects/Patek_project/force_data';


%% Simulation parameters

% Precision of simulation (default = 10^-5, use 10^-7 for high precision)
p.maxError = 10^-7;

% Duration of simulation (s)
p.simDur    = 0.010;

% Time values to evaluate results(use 1000-5000)
p.t = linspace(0,p.simDur,1000);


%% General parameter values
% Fixed for all individuals

% water I (kg m^2)
%p.waterI      = 1e-8;

% 3rd moment of area of dactyl (m^5) (for 12 mm long dactyl, r = 3 mm)
%p.dacA         = 0.25*(3e-3)*(12e-3)^4;

% Density of water (kg/m^3)
p.rho = 998;

%TODO: Adjust these to be custom for each individual


%% Run simulation


% Individual-specific parameter values
p = get_params(indiv,p);

% Knock out drag
%p.D = 0;

p.L3 = p.L3 +.1e-3;

%p.D = 3*p.D;

% Run simulation
[d,result] = run_sim(p);



%% Visualize simulation results
figure;

% subplot(3,1,1)
% plot(d.t.*1000,d.thetaIn.*180/pi,'r',d.t.*1000,d.thetaOut*180/pi,'b')
% legend('in ','out','Location','NorthWest');
% ylabel('angle (deg)')
% xlabel('time (ms)')
% grid on
% 
% subplot(3,1,2)
% plot(d.t.*1000,(d.dacAngSpd.*180/pi)./1000)
% ylabel('SB angular spd (deg/ms)')
% xlabel('time (ms)')
% grid on

subplot(2,1,1)
E_tot = d.E_kin + d.E_drag + d.E_elastic;
plot(d.t.*1000,d.E_elastic.*1000,'r',d.t.*1000,d.E_drag.*1000,'b',...
    d.t.*1000,d.E_kin.*1000,'g',d.t.*1000,E_tot.*1000,'k--')
ylabel('Energy (mJ)')
xlabel('time (ms)')
grid on
legend('elastic','drag','kinetic','total','Location','West');


subplot(2,1,2)
plot(d.t(2:end).*1000,diff(d.E_elastic)./diff(d.t),'r',...
    d.t(2:end).*1000,diff(d.E_drag)./diff(d.t),'b',...
    d.t(2:end).*1000,diff(d.E_kin)./diff(d.t),'g')
ylabel('Power (J/s)')
xlabel('time (ms)')
grid on


clear E_tot



%TODO: calc amd compare momentum values
%TODO: Set up to run all
%TODO: 



