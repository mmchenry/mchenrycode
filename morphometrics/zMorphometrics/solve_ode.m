function solve_ode
% Dynamic modeling of prey capture


%% Paths to data

zBaseM = '/Volumes/workgroup/Manuscripts/Relative flow/zMorphometrics_data';
dPath = '/Volumes/Docs/Projects/Relative velocity/lit_data';


%% Parameter values

% Density of water (kg m^-3)
rho_water   = 998; 
p.rho_water = rho_water;

%% Load data

% Load pooled data, mP
load([zBaseM filesep 'body_metrics'])

% Load data for bluegill, bass
load([dPath filesep 'lit_data_part2.mat'])

% Add gape diameter data
bass.gape = 2.93e-2;
blue.gape = 1.12e-2;


%% Define general parameters from loaded data

% Identify those w/out a SB
idx = mP.sb_vol==0;

% Define mean body density (no SB)
rho_noSB = mean(mP.rho_body(idx));

% Define mean body density (no SB)
rho_wSB = mean(mP.rho_body(~idx));

% Define I (with SB)
I_wSB = mean(mP.I(~idx));

% Define body volume (no SB)
V_noSB = mean(mP.Vbody(idx));

% Define body volume (with SB)
V_wSB = mean(mP.Vbody(~idx));

% Define mass (with SB)
M_wSB = mean(mP.Mbody(~idx));

% Mean COM with SB
COM_wSB = mean(mP.COM(~idx,:));

% Distance fromCOM and tail tip
tipLen = mean(mP.b_length(~idx))-mean(mP.COM(~idx,2));

% Mean position of level arm btn COV & COM with swim bladder
L = [mean(mP.COV(~idx,1)-mP.COM(~idx,1)) ...
    mean(mP.COV(~idx,2)-mP.COM(~idx,2)) ...
    mean(mP.COV(~idx,1)-mP.COM(~idx,1))];

clear mP



%% Sim 1: Bass

c2 = [0.098 -0.7 1.86 -2.19 1]; 

% Load up input parameters 
p.gapeDiameter  = bass.gape;
p.xInit         = 2.5e-3%2*p.gapeDiameter;
p.sp            = bass.cs;
p.rho_body      = rho_water;
p.tMax          = max(bass.t);
p.Uinit         = fnval(p.sp,0) .* polyval(c2,p.xInit.*1000);

% Run simulation
[T,X] = ode45(@(t,x) state_eqn(t,x,p),[0 p.tMax],[p.xInit p.Uinit 0]);

% Plot results
plotResults(T,X,p)


return


%% Sim 2: Bluegill


% Load up input parameters 
p.gapeDiameter  = blue.gape;
p.xInitial      = 10*p.gapeDiameter;
p.sp            = blue.cs;
p.rho_body      = rho_water;
p.tMax          = 0.9*max(blue.t);


% Run simulation
[T,X] = ode45(@(t,x) state_eqn(t,x,p),[0 p.tMax],[p.xInitial p.Uinit 0]);

% Plot results
plotResults(T,X,p)






function plotResults(T,X,p)
% Plots results of prey capture simulation



% From Day et al (2005)
%c2 = [0.348 -2.49 6.61 -7.78 3.56];

% Cut out data after prey enters the mouth
idx = X(:,1)>0;
t = T(idx);
X = X(idx,:);

% Translate results from X
bodyPos = X(:,1);
bodyVel = X(:,2);
relVel  = diff(X(:,3));

gapeDiameter = p.gapeDiameter;

clear T X idx


% Flow speed at the mouth over time
U_mouth = fnval(p.sp,t);

% Polynominal coefficients from Wainwright and Day (2007):
% flow velocity as a function of distance from mouth
c2 = [0.098 -0.7 1.86 -2.19 1]; 

% Flow speed at body position (m/s)
U_local = U_mouth.*polyval(c2,(bodyPos.*1000));


%local_velocity(U_mouth,x,gapeDiameter)

% Plot
figure
subplot(4,1,1)
plot(t,U_mouth); 
ylabel('Velocity at mouth (m/s)')
grid on

subplot(4,1,2)
plot(t,U_local);
hold on
plot(t,bodyVel,'r--')
ylabel('Local velocity (m/s)')
grid on
hold off

subplot(4,1,3)
plot(t,bodyPos,'-')
ylabel('x (mm)')
grid on

subplot(4,1,4)
plot(t(1:end-1),relVel,'-')
grid on
ylabel('Relative velocity (mm/s)')


function xdot = state_eqn(t,x,p)
% Differential equation defining the kinematics of a prey during suction
% feeding

% Translate inputs
bodyPos = x(1);
bodyVel = x(2);
relPos  = x(3);

% Body position relative to gape diameter
%bodyPos_rel = bodyPos;%/p.gapeDiameter;

% Spline coefficients
sp = p.sp;

% Water density 
rho_water = p.rho_water;

% Body density
rho_body = p.rho_body;

% From Day et al (2005)
%c2 = [0.348 -2.49 6.61 -7.78 3.56];

gapeDiameter = p.gapeDiameter;

clear x p


% Flow speed at the mouth at current instant (m/s)
U_mouth = fnval(sp,t);

% Flow acceleration at the mouth at current instant (m/s^2)
a_mouth = fnval(fnder(sp),t);

% Polynominal coefficients from Wainwright and Day (2007):
% flow velocity as a function of distance from mouth
cU = U_mouth.*[0.098 -0.7 1.86 -2.19 1]; 

cA = a_mouth.*[0.098 -0.7 1.86 -2.19 1]; 
    
% Flow speed at body position (m/s)
U_local = polyval(cU,bodyPos.*1000);

% Acceleration at body position (m/s^2)
a_local =  polyval(cA,bodyPos.*1000);

% Spatial rate of change in U at body position (1/s)
dUdx = polyval(polyder(cU),bodyPos.*1000);

% Local pressure gradient
dpdx  = rho_water .* (a_local + U_local.*dUdx);
%dpdx  = rho_water .* (a_local);

% Body velocity
xdot(1,1)  = bodyVel;

% Body acceleration
xdot(2,1)  = dpdx/rho_body;

% Relative acceleration
xdot(3,1) = dpdx/rho_water - dpdx/rho_body;




function U_local = local_velocity(U_mouth,x,gapeDiameter)
% From eqn. 34 in Muller

h1 = gapeDiameter/2;

U_local = U_mouth .* h1.^3 ./ ((x.^2 + h1.^2).^(3/2));


function dUdx = vel_gradient(U_mouth,x,gapeDiameter)
% First derivative wrt x of eqn. 34 from Muller

h1 = gapeDiameter/2;

dUdx = -(3.*U_mouth.*h1.^3.*x)./(h1.^2 + x.^2).^(5/2);



