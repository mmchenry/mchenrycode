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

dL = 1e-4;

indiv = 120;



%% General parameter values
% Fixed for all individuals

%TODO: Adjust these to be custom for each individual

% water I (kg m^2)
%p.waterI      = 1e-8;

% 3rd moment of area of dactyl (m^5) (for 12 mm long dactyl, r = 3 mm)
%p.dacA         = 0.25*(3e-3)*(12e-3)^4;

% Density of water (kg/m^3)
p.rho = 998;



%% Run for L3 

p = get_params(indiv,p);
p_start = p;


L3s = linspace(p.L3-dL,p.L3+dL,5);

clrs = {'r-','g-','b-'};



for i = 1:length(L3s)
    
    % Current L3 length
    p.L3 = L3s(i);
    
    %p.D = 2*p.D;
    
    %p.dacMass = p.dacMass*1.5;
    
    % Check geometry
    L = check_linkage(p,0);
    
    % Run the model with drag
    dD = run_sim(p,echoData);
    
    % Find index where max velocity is reached
    idx = find(dD.dacAngSpd == max(dD.dacAngSpd));
    
    % Store results
    rD.L3(i)         = p.L3;
    rD.alpha(i)      = max(dD.dacAngSpd);
    rD.E_drag(i)     = max(dD.E_drag);
    rD.theta_max(i)  = dD.thetaOut(idx)-dD.thetaOut(1);
    
    
    % Run the model without drag
    p.D = 0;
    dN = run_sim(p,echoData);
    
    % Find index where max velocity is reached
    idx = find(dD.dacAngSpd == max(dD.dacAngSpd));
    
    % Store results
    rN.L3(i)         = p.L3;
    rN.alpha(i)      = max(dN.dacAngSpd);
    
    rN.theta_max(i)  = dN.thetaOut(idx)-dN.thetaOut(1);

%     r.al_pred(i)     = sqrt(p.kSpring./((p.dacI+p.waterI))) * ...
%                              ((p.thetaRest-p.thetaStart));
    
    % Reset parameter values for next loop
    p = p_start;

end

figure;
subplot(3,1,1)
plot(rD.L3,(180/pi).*rD.theta_max,'o-',...
     rN.L3,(180/pi).*rN.theta_max,'ro-')
ylabel('Range of motion (deg)')
xlabel('L3')
legend('w/drag','no drag')
axis square

subplot(3,1,2)
plot(rD.L3,1000.*rD.E_drag,'o-')
ylabel('Energy loss (mJ)')
xlabel('L3')
axis square

subplot(3,1,3)
plot(rD.L3,rD.alpha.*(180/pi),'o-',...
     rN.L3,rN.alpha.*(180/pi),'ro-')
ylabel('max ang speed (deg/s)')
xlabel('L3')
axis square

return

figure
subplot(3,1,1)
plot(r.L3,r.al_nodrag,'o',r.L3,r.al_pred,'-',r.L3,r.al_drag,'or')
subplot(3,1,2)
plot(r.L3,r.spring_disp,'o-')
subplot(3,1,3)
plot(r.L3,r.max_kine,'o-')


ttt=1