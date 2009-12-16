function [d,result] = run_sim(p,echoData)
% This collects parameters values, saves them to disk for the Mathematica
% kernel to find, runs the kernel, reads the resutls from disk, and passes
% them to the output structure d.  Takes the following inputs:
%
% txtPath    - (string) Path to text file of Mathematica commands
% kernelPath - (string) Path to the Mathematica kernel
% currPath - (string)   Current path (stores info on current sim)
%

result = {};

if nargin < 2
    echoData = 1;
end

%% Check input geometry - return empty "d", if geometry not possible

L = check_linkage(p,0);

% Report problems with the geometry
if (p.thetaStart < L.thetaInMin)
    result{length(result)+1} = ...
     'Starting thetaIn value below what is possible for linkage geometry';

elseif (p.thetaStart > L.thetaInMax)
   result{length(result)+1} = ...
     'Starting thetaIn value above what is possible for linkage geometry';  
end

if isempty(L)
    result{length(result)+1} = ...
     'No range of motion possible for requested geomtry';  
end

% Stop execution if there is a reported problem
if ~isempty(result)
    d = [];
    warning(result)
    disp('No simulation run')
    return
end

clear h_BD si h_BF gamma A B C D E F G 


%% Export model parameters for mathematica   

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
save([p.simsPath filesep 'input_params.mat'],'data','-v4')

% Save time vector for evaluating results of model
time = p.t;
save([p.simsPath filesep 'eval_time.mat'],'time','-v4')

clear time data


%% Run sim code with mathematica kernel

if echoData
    disp(' ')
    disp(['Executing ' txtFile ' . . .'])
    tic
end

[status,res_text] = unix([p.kernelPath ' -noprompt -run ' p.mathFile]);

if echoData
    disp(['. . . done (' num2str(toc) ' s)']);
end

%% Load and return simulation results

% Load coordinate data (nx6 matrix of 2d values for pos., vel., accel.)
carp1   = importMathematica([p.simsPath filesep 'carpusP1.mat']);
carp2   = importMathematica([p.simsPath filesep 'carpusP2.mat']);
dac     = importMathematica([p.simsPath filesep 'dactylP1.mat']);
mV      = importMathematica([p.simsPath filesep 'mVP1.mat']);
gnd     = importMathematica([p.simsPath filesep 'groundP1.mat']);

% Load angular data
%mV_ang  = importMathematica([p.simsPath filesep 'mvAng.mat']);
%dac_ang = importMathematica([p.simsPath filesep 'dacAng.mat']);

% Load torque data
springTau = importMathematica([p.simsPath filesep 'springMoment.mat']);
dragTau   = importMathematica([p.simsPath filesep 'dragMoment.mat']);

% Calculate mV angular position
mVAng  = unwrap(atan2(mV(:,2),mV(:,1)));

% Calculate input angle
thetaIn  = unwrap(pi/2 - mVAng);

% Calculate output angle
h_BD  = sqrt(p.L1^2 + p.L2^2 - 2*p.L1*p.L2*cos(thetaIn));
thetaOut = real(acos((p.L3^2+p.L4^2-h_BD.^2)./(2*p.L3*p.L4)));

dacAng = unwrap(atan2(dac(:,2)-carp2(:,2),dac(:,1)-carp2(:,1)));

% Trim data beyond where output angle approaches maximum range of motion
if max(thetaOut > 0.97*L.thetaOutMax) > 0
    idx = 1:find(thetaOut > 0.97*L.thetaOutMax,1);
else
    idx = 1:length(thetaOut);
end

% Define time vectors
t_d = [p.t(idx)]';
t_v = t_d(1:end-1,1);
t_a = t_d(2:end-1,1);

% Store time vector for results
d.t = t_a;

% Calculate angular velocity
dacAngSpd = diff(dacAng(idx))./diff(t_d);

% Calculate dactyl speed
dacCOMSpd      = sqrt(diff(dac(idx,1)).^2 + diff(dac(idx,2)).^2)./diff(t_d);
dacBaseSpd     = sqrt(diff(carp2(idx,1)).^2 + diff(carp2(idx,2)).^2)./diff(t_d);

%d.mVAngVel      = diff(d.mVAng)./diff(t_d);
%d.mVAngAccel    = diff(d.mVAngVel)./diff(t_v);

% Interpolate angular data
d.thetaIn   = interp1(t_d,thetaIn(idx),t_a);
d.thetaOut  = interp1(t_d,thetaOut(idx),t_a);

% Interpolate speed data
d.dacAngSpd  = interp1(t_v,dacAngSpd(idx(1:end-1)),t_a);
d.dacCOMSpd  = interp1(t_v,dacCOMSpd(idx(1:end-1)),t_a);
d.dacBaseSpd = interp1(t_v,dacBaseSpd(idx(1:end-1)),t_a);

% Interpolate torque data
d.springTau   = interp1(t_d,springTau(idx),t_a);
d.dragTau     = interp1(t_d,dragTau(idx),t_a);

% Calculate elastic energy (assumes negative displacement of mV)
d.E_elastic  = abs(0.5 .* d.springTau .* (d.thetaIn-p.thetaRest));

% Calculate drag energy
d.E_drag = cumtrapz(abs(d.thetaOut(1)-d.thetaOut),abs(d.dragTau));

% Calculate kinetic energy
d.E_kin = (0.5 * (p.dacI+p.waterI) .* d.dacAngSpd.^2) + ...
          (0.5 * p.dacMass .* d.dacBaseSpd.^2);
    
% Interpolate position/velocity/acceleration data
d.carp1PVA  = interp1(t_d,carp1(idx,:),t_a);
d.carp2PVA  = interp1(t_d,carp2(idx,:),t_a);
d.dacPVA    = interp1(t_d,dac(idx,:),t_a);
d.mVPVA     = interp1(t_d,mV(idx,:),t_a);
d.gndPVA    = gnd;
d.KT        = (p.L1*p.L2)/(p.L3*p.L4) .* csc(d.thetaOut) .* ...
              sqrt(1-((p.L1^2+p.L2^2-p.L3^2-p.L4^2+2*p.L3*p.L4.*cos(d.thetaOut)).^2) ...
              ./ (4*p.L1^2 *p.L2^2));


%% Delete temporary files from disk

delete([p.simsPath filesep 'input_params.mat'])
delete([p.simsPath filesep 'eval_time.mat'])
delete([p.simsPath filesep 'carpusP1.mat'])
delete([p.simsPath filesep 'carpusP2.mat'])
delete([p.simsPath filesep 'dactylP1.mat']);
delete([p.simsPath filesep 'mVP1.mat']);
delete([p.simsPath filesep 'springMoment.mat']);
delete([p.simsPath filesep 'dragMoment.mat']);
delete([p.simsPath filesep 'groundP1.mat']);


%% Update result text
result{1} = res_text;


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



