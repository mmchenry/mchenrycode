function d = run_sim(p,txtFile,kernelPath,simsPath)
% This collects parameters values, saves them to disk for the Mathematica
% kernel to find, runs the kernel, reads the resutls from disk, and passes
% them to the output structure d.  Takes the following inputs:
%
% txtPath    - (string) Path to text file of Mathematica commands
% kernelPath - (string) Path to the Mathematica kernel
% currPath - (string)   Current path (stores info on current sim)
%

%% Export model parameters for mathematica   

% Check input geometry - return empty "d", if geometry not possible
h  = sqrt(p.L1^2 + p.L2^2 - 2*p.L1*p.L2*cos(p.thetaStart));
%     si = acos((h^2+p.L1^2-p.L2^2)/(2*h*p.L1)) + ...
%          acos((h^2+p.L4^2-p.L3^2)/(2*h*p.L4));

if (p.L4 > (p.L3 + h)) || (p.L4 < h-p.L3)
    d = [];
    return
end

% The order of these inputs must be consistent with 'sim_code.txt'
data(1)  = p.L1;
data(2)  = p.L2;
data(3)  = p.L3;
data(4)  = p.L4;
data(5)  = p.L5;
data(6)  = p.h_AF;
data(7)  = p.L_COM;
data(8)  = 0;
data(9)  = p.thetaStart;
data(10) = p.dacMass;
data(11) = p.dacI;
data(12) = p.kSpring;
data(13) = p.thetaRest;
data(14) = p.simDur;
data(15) = p.dacA;
data(16) = p.rho;
data(17) = p.Cd;
data(18) = p.waterI;
data(19) = p.maxError;

% Save data for use by mathematica model
save([simsPath filesep 'input_params.mat'],'data','-v4')

% Save time vector for evaluating results of model
time = p.t;
save([simsPath filesep 'eval_time.mat'],'time','-v4')

clear time data


%% Run sim code

disp(' ')
disp(['Executing ' txtFile ' . . .'])
tic
[status,result] = unix([kernelPath ' -noprompt -run ' txtFile]);
disp(result)
disp(['. . . done (' num2str(toc) ' s)']);


%% Load and return simulation results

% Load postion data (nx6 matrix of 2d values for pos., vel., accel.)
carp1   = importMathematica([simsPath filesep 'carpusP1.mat']);
carp2   = importMathematica([simsPath filesep 'carpusP2.mat']);
dac     = importMathematica([simsPath filesep 'dactylP1.mat']);
mV      = importMathematica([simsPath filesep 'mVP1.mat']);

% Load torque data
d.springTau = importMathematica([simsPath filesep 'springMoment.mat']);
d.dragTau   = importMathematica([simsPath filesep 'dragMoment.mat']);

% Define time vectors
t_d = [p.t]';
t_v = t_d(1:end-1,1);
t_a = t_d(2:end-1,1);

% Store time vector to evaluate results
d.t = t_a;

% Calculate mV angular displacement, vel., accel.
d.mVAng         = atan2(mV(:,2),mV(:,1));
d.mVAngVel      = diff(d.mVAng)./diff(t_d);
d.mVAngAccel    = diff(d.mVAngVel)./diff(t_v);

% Interpolate mV data to t_a
d.mVAng     = interp1(t_d,d.mVAng,t_a);
d.mVAngVel  = interp1(t_v,d.mVAngVel,t_a);

% Calculate dactyl kinematics
d.dacAng      = atan2(dac(:,2)-carp2(:,2),dac(:,1)-carp2(:,1));
d.dacAngVel   = diff(d.dacAng)./diff(t_d);
d.dacAngAccel = diff(d.dacAngVel)./diff(t_v);
d.dacSpd      = sqrt(diff(carp2(:,1)).^2 + diff(carp2(:,2)).^2)./diff(t_d);

% Interpolate dactyl data wrt t_a
d.dacAng     = interp1(t_d,d.dacAng,t_a);
d.dacAngVel  = interp1(t_v,d.dacAngVel,t_a);
d.dacSpd     = interp1(t_v,d.dacSpd,t_a);

% Interpolate torque data
d.springTau   = interp1(t_d,d.springTau,t_a);
d.dragTau     = interp1(t_d,d.dragTau,t_a);

% Calculate elastic energy (assumes negative displacement of mV)
d.E_elastic  = abs(0.5 .* d.springTau .* (pi/2-d.mVAng-p.thetaRest));

% Calculate drag energy
d.E_drag = cumtrapz(abs(d.dacAng(1)-d.dacAng),abs(d.dragTau));

% Calculate kinetic energy
d.E_kin = (0.5 * (p.dacI+p.waterI) .* d.dacAngVel.^2) + (0.5 * p.dacMass .* d.dacSpd.^2);
    
% Interpolate position/velocity/acceleration data
d.carp1PVA  = interp1(t_d,carp1,t_a);
d.carp2PVA  = interp1(t_d,carp2,t_a);
d.dacPVA    = interp1(t_d,dac,t_a);
d.mVPVA     = interp1(t_d,mV,t_a);


%% Delete temporary files

delete([simsPath filesep 'input_params.mat'])
delete([simsPath filesep 'eval_time.mat'])
delete([simsPath filesep 'carpusP1.mat'])
delete([simsPath filesep 'carpusP2.mat'])
delete([simsPath filesep 'dactylP1.mat']);
delete([simsPath filesep 'mVP1.mat']);
delete([simsPath filesep 'springMoment.mat']);
delete([simsPath filesep 'dragMoment.mat']);


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



