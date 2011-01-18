function vary_linkage


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
p.simDur    = 0.015;

% Time values to evaluate results(use 1000-5000)
p.t = linspace(0,p.simDur,1000);


visData = 0;
echoData = 0;
animate = 0;

dL = .1e-3;

indiv = 120;



%% General parameter values
% Fixed for all individuals

%TODO: Adjust these to be custom for each individual

% water I (kg m^2)
p.waterI      = 1e-8;

% 3rd moment of area of dactyl (m^5) (for 12 mm long dactyl, r = 3 mm)
p.dacA         = 0.25*(3e-3)*(12e-3)^4;

% Density of water (kg/m^3)
p.rho = 998;



%% Run for L3 

p = get_params(indiv,p);
p_start = p;


L3s = [p.L3-dL p.L3 p.L3+dL];

clrs = {'r-','g-','b-'};



for i = 1:length(L3s)
    
    % Current L3 length
    p.L3 = L3s(i);
    
    % Check geometry
    L = check_linkage(p,0);
    
    % Run the model with drag
    %dD = run_sim(p,echoData);
    
    % Run the model without drag
    p.D = 0;
    dN = run_sim(p,echoData);
    
    r.L3(i)          = p.L3;
    r.al_nodrag(i)   = max(dN.dacAngSpd);
    r.al_pred(i)     = sqrt(p.kSpring./((p.dacI+p.waterI) + p.dacMass*p.L5^2)) * ...
                             (p.thetaRest-p.thetaStart);
    
    % Reset parameter values for next loop
    p = p_start;

end

figure
plot(r.L3,r.al_nodrag,'o',r.L3,r.al_pred,'-')


ttt=1
