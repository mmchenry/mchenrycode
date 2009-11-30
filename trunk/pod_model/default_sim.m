function default_sim
% Runs a simulation with the default parameter values for falcatus


%% Components of code to execute after simulation
visInput    = 0; % Plots input coordinates
animate     = 1; % Animates the simulation results
visTorques  = 0; % Examines the forces predicted by the model
visKine     = 0;
visEnergy   = 0;

%% Define Paths

% Path to text file of Mathematica commands
txtFile = ['"<</Volumes/Docs/Projects/Patek_project/pod_model/sim_code.txt"'];

% Path to the Mathematica kernel
kernelPath = '/Applications/Mathematica.app/Contents/MacOS/MathKernel';

% Root path for simulation data
simsPath = '/Volumes/Docs/Projects/Patek_project/sims';

% Define and save current sim path
currPath = [simsPath filesep '/s001'];

% Save current path in text file
%dlmwrite([simsPath filesep 'currPath.txt'],currPath,'delimiter','%s')


%% Define model parameter set 

% Linkage coordinates (m)
p.L1        = 4.5e-3;
p.L2        = 3.0e-3; 
p.L3        = 1.0e-3;
p.L4        = 3.8e-3;
p.L5        = 1.5e-3;
%p.L6        = 4.0e-3;
p.h_AF        = 0.5e-3;

% Distance from E to the COM
p.L_COM     = 5.0e-3;

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

%p.currPath = currPath;

% Clear results from prior simulation within current path
%delete([currPath filesep '*.mat'])

% Run the model
d = run_sim(p,txtFile,kernelPath,simsPath);

if isempty(d)
    error('Impossible geometry')
end


%% Plot of coorindates describing the initial morphology

% Check that initial geometry matches the requested geometry 
% (The model has the ability to alter the geometry)


if visInput
    

    % Distance between B & D
    h_BD  = sqrt(p.L1^2 + p.L2^2 - 2*p.L1*p.L2*cos(p.thetaStart));
    
    % Angl btwn L4 & L1
    si = acos((h_BD.^2+p.L1^2-p.L2^2)/(2*h_BD*p.L1)) + ...
         acos((h_BD^2+p.L4^2-p.L3^2)/(2*h_BD*p.L4));
    
    % Distance between points B & F
    h_BF    = sqrt(p.h_AF^2 + p.L2^2 + 2*p.h_AF*p.L2*cos(p.thetaStart));
    
    % Angle btwn 5 and the x-axis
    gamma = pi - acot(p.L5/sqrt(h_BF^2-p.L5^2))...
               - atan(cot(p.thetaStart) + p.h_AF*csc(p.thetaStart)/p.L2);
    
       
    % Define points
    A = [0 0];
    B = [p.L2*sin(p.thetaStart) p.L2*cos(p.thetaStart)];
    C = [p.L4*sin(si) p.L1-p.L4*cos(si)];
    D = [0 p.L1];
    E = [B(1)+p.L5*cos(gamma) ...
         B(2)-p.L5*sin(gamma)];
    F = [0 -p.h_AF];
    G = [E(1) - p.L_COM * sin(gamma)
         E(2) - p.L_COM * cos(gamma)];
       
    figure;
    
    plot([A(1) B(1)],[A(2) B(2)],'r',...
         [B(1) C(1)],[B(2) C(2)],'g',...
         [C(1) D(1)],[C(2) D(2)],'b',...
         [B(1) E(1)],[B(2) E(2)],'m',...
         [E(1) F(1)],[E(2) F(2)],'m-',...
         [E(1) G(1)],[E(2) G(2)],'k--')

     axis equal
   
    
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