function run_sim

%% Components of code to execute
updateParam     = 1; % Saves a new set of parameter files 
runSim          = 1; % Uses mathematica to run a simulation
animate         = 1; % Animates the simulation results



%% Define Paths

% Path to text file of Mathematica commands
txtFile = ['"<</Volumes/Docs/Projects/Patek_project/pod_model/sim_code.txt"'];

% Path to the Mathematica kernel
kerPath = '/Applications/Mathematica.app/Contents/MacOS/MathKernel';

% Root path for simulation data
simsPath = '/Volumes/Docs/Projects/Patek_project/pod_model/sims';

% Define and save current sim path
currPath = [simsPath filesep '/s001'];

% Save current path in text file
dlmwrite([simsPath filesep 'currPath.txt'],currPath,'delimiter','%s')

% Clear results from prior simulation within current path
delete([currPath filesep '*.mat'])


%% Define model parameter set 

if updateParam
    % Linkage coordinates (m)
    p.L1        = 5.26e-3;
    p.L2        = 3.18e-3;
    p.L3        = 0.81e-3;
    p.L4        = 5.11e-3;
    p.EXLocal   = 1.43e-3;
    p.EYLocal   = 1.64e-3;
    p.FXLocal   = 5.52e-3;
    p.FYLocal   = 4.30e-3;
    
    % Initial input angle (rad)
    p.theta0    = (66/180)*pi;
    
    % dactyl mass (kg)
    p.dMass     = 4.53e-4;
    
    % dactyl I (kg m^2)
    p.dI        = 8.16e-8;
    
    % Torsion spring at mV joint (Nm/rad)  (= linear spring stiffness (N/m) * L2^2
    % (approx. distance to force application)
    p.kSpring   = (60*10^3) * (3.18 *10^-3)^2;
    
    % Maximum torque from the spring (Nm), estimated as max force (29 N in Zack et
    % al paper), times L2 (3.18e-3, approximate distance to force application)
    p.Tmax      =  29 * (3.18 *10^-3);
    
    % Duration of simulation (s)
    p.simDur    = 0.01;
    
    % Time values to evaluate results
    p.t = linspace(0,p.simDur,100);
    
    save([currPath filesep 'param_struct.mat'],'p')
    
    clear p
end

%% Export model parameters for mathematica

if updateParam
    
    % Load p, the structure of parameter values
    load([currPath filesep 'param_struct.mat'])
    
    % The order of these inputs must be consistent with 'sim_code.txt'
    data(1)  = p.L1;
    data(2)  = p.L2;
    data(3)  = p.L3;
    data(4)  = p.L4;
    data(5)  = p.EXLocal;
    data(6)  = p.EYLocal;
    data(7)  = p.FXLocal;
    data(8)  = p.FYLocal;
    data(9)  = p.theta0;
    data(10) = p.dMass;
    data(11) = p.dI;
    data(12) = p.kSpring;
    data(13) = p.Tmax;
    data(14) = p.simDur;
    
    % Save
    save([currPath filesep 'input_params.mat'],'data','-v4')
    
    % Time vector for evaluating results
    time = p.t;
    save([currPath filesep 'eval_time.mat'],'time','-v4')
    
    clear time data p
end

%% Run sim code

if runSim
    disp(' ')
    disp(['Executing ' txtFile ' . . .'])
    tic
    [status,result] = unix([kerPath ' -noprompt -run ' txtFile]);
    disp(result)
    disp(['. . . done (' num2str(toc) ' s)']);
end


%% Import simulation results

% Postion data are a nx6 matrix of 2d values for position, velocity, accel

carp1   = importMathematica([currPath filesep 'carpusP1.mat']);
carp2   = importMathematica([currPath filesep 'carpusP2.mat']);
dac     = importMathematica([currPath filesep 'dactylP1.mat']);
mV      = importMathematica([currPath filesep 'mVP1.mat']);


%% Visualize position data 

% figure;
% plot(p.t,carp1(:,1),'-',p.t,carp1(:,2),'r-')


%% Animate simulation

if animate
    figure;
    
    for i = 1:length(p.t)
        h1 = plot([0 mV(i,1)],[0 mV(i,2)],'b-');
        hold on
        axis square
        ylim([-10e-3,3e-3])
        xlim([-3e-3,10e-3])
        h2 = plot([carp1(i,1) carp2(i,1)],[carp1(i,2) carp2(i,2)],'r-');
        h3 = plot([carp2(i,1) dac(i,1)],[carp2(i,2) dac(i,2)],'g-');
        title(['t = ' num2str(p.t(i))])
        pause(.01)
        hold off
    end
end


return

%save([sims_path filesep 'curr_path.mat'],'curr_path','-v4')

% Execute Mathematica sim code




%<</Volumes/Docs/Projects/Patek_project/Mathematica/iotest.txt

function data = importMathematica(filePath)

% Check that no "Expression" item is in the workspace
a = whos('Expression*');
if ~isempty(a)
    error('Clear Expression vector in the workspace before importing')
end
clear a

% Load Mathematica file & redefine new "Expression" matrix
load(filePath)
a = whos('Expression*');
if length(a) > 1
    error('More than one Expression vector in the workspace')
elseif isempty(a)
    error('No Expression vector in the workspace')
end

eval(['data = ' a.name ';']);
data = data';



