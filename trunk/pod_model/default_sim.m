function default_sim
% Runs a simulation with the default parameter values for falcatus


%% Components of code to execute after simulation
visInput    = 1; % Plots input coordinates
animate     = 0; % Animates the simulation results
visTorques  = 1; % Examines the forces predicted by the model
visKine     = 1;
visEnergy   = 1;

%% Define Paths

% Path to text file of Mathematica commands
txtFile = ['"<</Volumes/Docs/Projects/Patek_project/pod_model/sim_code.txt"'];

% Path to the Mathematica kernel
kernelPath = '/Applications/Mathematica.app/Contents/MacOS/MathKernel';

% Root path for simulation data
simsPath = '/Volumes/Docs/Projects/Patek_project/pod_model/sims';

% Define and save current sim path
currPath = [simsPath filesep '/s001'];

% Save current path in text file
dlmwrite([simsPath filesep 'currPath.txt'],currPath,'delimiter','%s')


%% Define model parameter set 

% Linkage coordinates (m)
p.L1        = 5.26e-3;
p.L2        = 2.2e-3; %3.18e-3;
p.L3        = 0.81e-3;
p.L4        = 4.8e-3;
p.EXLocal   = 1.43e-3;
p.EYLocal   = 1.64e-3;
p.FXLocal   = 5.52e-3;
p.FYLocal   = 4.30e-3;

% Initial input angle (rad)
p.thetaStart  = (55/180)*pi;

% dactyl mass (kg)
p.dacMass     = 4.53e-4;

% dactyl I (kg m^2)
p.dacI        = 8.16e-8;

% water I (kg m^2
p.waterI      = 1e-8;

% 3rd moment of area of dactyl (m^5) (for 12 mm long dactyl, r = 3 mm)
p.dacA         = 0.25*(3e-3)*(12e-3)^4;

% Drag coefficient of elliptical dactyl element (Approx)
p.Cd           = 0*1;

% Torsion spring at mV joint (Nm/rad)  (= linear spring stiffness (N/m) * L2^2
% (approx. distance to force application)
p.kSpring   = (60*10^3) * (3.18 *10^-3)^2;

% Resting position for torsion spring
p.thetaRest = (85/180)*pi;

% Density of water (kg/m^3)
p.rho = 998;

% Duration of simulation (s)
p.simDur    = 0.004;

% Time values to evaluate results(use 1000-5000)
p.t = linspace(0,p.simDur,1000);

% Precision of simulation (default = 10^-5, 10^-7 for high precision)
p.maxError = 10^-7;

p.currPath = currPath;

% Clear results from prior simulation within current path
delete([currPath filesep '*.mat'])

% Run the model
d = run_sim(p,txtFile,kernelPath,currPath,simsPath);

if isempty(d)
    error('Impossible geometry')
end


%% Plot of coorindates describing the initial morphology

% Check that initial geometry matches the requested geometry 
% (The model has the ability to alter the geometry)




if visInput
    
    
    hTmp = sqrt(p.L3^2 + p.L6^2)
    
    %h  = sqrt(p.L1^2 + p.L2^2 - 2*p.L1*p.L2*cos(p.thetaStart));
    %si = acos((h.^2+p.L1^2-p.L2^2)/(2*h*p.L1));
    
    
end

%% Animate simulation

if animate
    
    % Parameters
    repeat = 3;       % Number of runs
    maxFrames = 100;  % Max number of frames in animation
    
    skipFrame = round(size(d.dacPVA,1)/maxFrames);
    
    hf = figure;
    set(hf,'DoubleBuffer','on')
    pause(.01)
    
    for j = 1:repeat
        for i = 1:skipFrame:length(d.t)
            h1 = plot([0 d.mVPVA(i,1)],[0 d.mVPVA(i,2)],'k-');
            hold on
            axis square
            ylim([-11e-3,6e-3])
            xlim([-3e-3,14e-3])
            h2 = plot([d.carp1PVA(i,1) d.carp2PVA(i,1)],...
                      [d.carp1PVA(i,2) d.carp2PVA(i,2)],'r-');
            h3 = plot([d.carp2PVA(i,1) d.dacPVA(i,1)],...
                      [d.carp2PVA(i,2) d.dacPVA(i,2)],'g-');
            h4 = plot([d.mVPVA(i,1) d.carp1PVA(i,1)],...
                      [d.mVPVA(i,2) d.carp1PVA(i,2)],'b--');
            h5 = plot([0 d.carp1PVA(i,1)],...
                      [p.L1 d.carp1PVA(i,2)],'b--');
            h6 = plot([0 0],[0 p.L1],'b--');
            title(['t = ' num2str(d.t(i))])
            pause(.02)
            hold off
            
            clear h1 h2 h3 h4 h5 h6
        end
    end
    
    clear skipFrame maxFrames repeat hf i j
end


%% Analyze predicted force/moments

if visTorques
 
    figure;
    plot(d.t,d.springTau,'-')
    grid on
    hold on
    plot(d.t,d.dragTau,'r-')
    grid on
    hold off
    ylabel('Torque (N m)')
    xlabel('Time (s)')
    legend('spring','drag','Location','SouthWest')
end


%% Energetic analysis

if visEnergy
    
    figure;
    plot(d.t,d.E_elastic,'-')
    hold on
    plot(d.t,d.E_kin,'r-')
    plot(d.t,d.E_drag,'g-')
    plot(d.t,d.E_kin+d.E_elastic,'k--')
    ylabel('Energy (j)')
    xlabel('time (s)')
    legend('elastic','kine','drag','elast+kin','Location','NorthWest')
    grid on
    hold off
    
end
                 

%% Visualize kinematics

if visKine
    
    figure;
    subplot(3,1,3);
    plot(d.t,d.mVAng.*(180/pi),'-')
    hold on
    plot(d.t,d.dacAng.*(180/pi),'r-')
    grid on
    hold off
    ylabel('Position (deg)')
    
    subplot(3,1,2)
    plot(d.t,d.mVAngVel.*(180/pi),'-')
    hold on
    plot(d.t,d.dacAngVel.*(180/pi),'r-')
    grid on
    hold off
    ylabel('Velocity (deg/s)')
    
    subplot(3,1,1)
    plot(d.t,d.mVAngAccel.*(180/pi))
    hold on
    plot(d.t,d.dacAngAccel.*(180/pi),'r-')
    grid on
    hold off
    ylabel('Acceleration (deg/s^2)')
    legend('mV','dactyl','Location','SouthWest')
    
end