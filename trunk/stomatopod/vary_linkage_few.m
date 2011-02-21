function vary_linkage_few
% Similar to vary_linkage, but run on just a few sequences


%% Simulation parameters

% Precision of simulation (default = 10^-5, use 10^-7 for high precision)
p.maxError = 10^-7;

% Duration of simulation (s)
p.simDur    = 0.02;

% Time values to evaluate results(use 1000-5000)
p.t = linspace(0,p.simDur,1000);


visData = 0;
echoData = 0;
animate = 0;
indiv = 122;

p = get_params(indiv,p);

L3_min = p.L3 - .0001;
L3_max = p.L3 + .0003;



%% General parameter values
% Fixed for all individuals

%TODO: Adjust these to be custom for each individual

% water I (kg m^2)
%p.waterI      = 1e-8;

% 3rd moment of area of dactyl (m^5) (for 12 mm long dactyl, r = 3 mm)
%p.dacA         = 0.25*(3e-3)*(12e-3)^4;

% Density of water (kg/m^3)
p.rho = 998;



clrs = {'r-','g-','b-'};

mksize = 4;

%% With drag, low L3

% Current L3 length
p = get_params(indiv,p);
p.L3 = L3_min;

% Check geometry
L = check_linkage(p,0);

% Run the model with drag
dD = run_sim(p,echoData);

% Find index where max velocity is reached
idx = find(dD.Dgamma == max(dD.Dgamma));

% Store results
lD.L3         = p.L3;
lD.theta      = dD.theta;
lD.gamma      = dD.gamma;
lD.Dgamma     = dD.Dgamma;
lD.E_drag     = dD.E_drag;
lD.theta_max  = dD.gamma(idx)-dD.gamma(1);
lD.t_peak     = dD.t(dD.Dgamma==max(dD.Dgamma));
lD.KT_min     = min(dD.KT);
lD.efficiency = 100*max(dD.E_kin)./max(dD.E_elastic);
 

%% With drag, high L3

% Current L3 length
p = get_params(indiv,p);
p.L3 = L3_max;

% Check geometry
L = check_linkage(p,0);

% Run the model with drag
dD = run_sim(p,echoData);

% Find index where max velocity is reached
idx = find(dD.Dgamma == max(dD.Dgamma));

% Store results
hD.L3         = p.L3;
hD.theta      = dD.theta;
hD.gamma      = dD.gamma;
hD.Dgamma     = dD.Dgamma;
hD.E_drag     = dD.E_drag;
hD.theta_max  = dD.gamma(idx)-dD.gamma(1);
hD.t_peak     = dD.t(dD.Dgamma==max(dD.Dgamma));
hD.KT_min     = min(dD.KT);
hD.efficiency = 100*max(dD.E_kin)./max(dD.E_elastic);


%% Without drag, low L3

% Current L3 length
p = get_params(indiv,p);
p.L3 = L3_min;
p.D = 0;

% Check geometry
L = check_linkage(p,0);

% Run the model with drag
dD = run_sim(p,echoData);

% Find index where max velocity is reached
idx = find(dD.Dgamma == max(dD.Dgamma));

% Store results
lN.L3         = p.L3;
lN.theta      = dD.theta;
lN.gamma      = dD.gamma;
lN.Dgamma     = dD.Dgamma;
lN.E_drag     = dD.E_drag;
lN.theta_max  = dD.gamma(idx)-dD.gamma(1);
lN.t_peak     = dD.t(dD.Dgamma==max(dD.Dgamma));
lN.KT_min     = min(dD.KT);
lN.efficiency = 100*max(dD.E_kin)./max(dD.E_elastic);
 

%% Without drag, high L3

% Current L3 length
p = get_params(indiv,p);
p.L3 = L3_max;
p.D = 0;

% Check geometry
L = check_linkage(p,0);

% Run the model with drag
dD = run_sim(p,echoData);

% Find index where max velocity is reached
idx = find(dD.Dgamma == max(dD.Dgamma));

% Store results
hN.L3         = p.L3;
hN.theta      = dD.theta;
hN.gamma      = dD.gamma;
hN.Dgamma     = dD.Dgamma;
hN.E_drag     = dD.E_drag;
hN.theta_max  = dD.gamma(idx)-dD.gamma(1);
hN.t_peak     = dD.t(dD.Dgamma==max(dD.Dgamma));
hN.KT_min     = min(dD.KT);
hN.efficiency = 100*max(dD.E_kin)./max(dD.E_elastic);
    


%% Plot results

figure;

subplot(1,2,1)
plot(lD.gamma,lD.Dgamma.^2,'b-',hD.gamma,hD.Dgamma.^2,'b--')
legend(num2str(lD.KT_min),num2str(hD.KT_min),'Location','NorthWest')
xlabel('gamma');
ylabel('Dgamma^2');
axis square


return

subplot(2,1,1)
h = plot(1000.*rD.L3,1000.*rD.t_peak,'o-',...
         1000.*rN.L3,1000.*rN.t_peak,'ro-');
legend('w/drag','no drag')
ylabel('time to peak speed (ms)')
xlabel('L3')
axis square

set(h(1),'MarkerSize',mksize)
set(h(1),'MarkerFaceColor','b')
set(h(2),'MarkerSize',mksize)
set(h(2),'MarkerFaceColor','r')

subplot(2,1,2)
h = plot(1000.*rD.L3,rD.gamma.*(180/pi)./1000,'o-',...
         1000.*rN.L3,rN.gamma.*(180/pi)./1000,'ro-',...
         1000.*rN.L3,rN.al_pred.*(180/pi)./1000,'k-');
     
set(h(1),'MarkerSize',mksize)
set(h(1),'MarkerFaceColor','b')
set(h(2),'MarkerSize',mksize)
set(h(2),'MarkerFaceColor','r')

ylabel('max ang speed (deg/ms)')
xlabel('L3')
axis square

figure
subplot(2,3,[1:2,4:5])
h = plot(rD.KT_min,rD.gamma.*(180/pi)./1000,'o-',...
         rN.KT_min,rN.gamma.*(180/pi)./1000,'ro',...
         rN.KT_min,rN.al_pred.*(180/pi)./1000,'k-');
     
set(h(1),'MarkerSize',mksize)
set(h(1),'MarkerFaceColor','b')
set(h(2),'MarkerSize',mksize)
set(h(2),'MarkerFaceColor','r')
ylabel('max ang speed (deg/ms)')
xlabel('KT min')

figure
subplot(2,3,[1:2,4:5])
h = plot(rD.KT_min,rD.efficiency,'-',...
         rN.KT_min,rN.efficiency,'r-');
     
set(h(1),'MarkerSize',mksize)
set(h(1),'MarkerFaceColor','b')
set(h(2),'MarkerSize',mksize)
set(h(2),'MarkerFaceColor','r')
ylabel('Efficiency')
xlabel('KT min')
%axis square


return

% Plot results
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
plot(rD.L3,rD.gamma.*(180/pi),'o-',...
     rN.L3,rN.gamma.*(180/pi),'ro-')
ylabel('max ang speed (deg/s)')
xlabel('L3')
axis square





figure
subplot(3,1,1)
plot(r.L3,r.al_nodrag,'o',r.L3,r.al_pred,'-',r.L3,r.al_drag,'or')
subplot(3,1,2)
plot(r.L3,r.spring_disp,'o-')
subplot(3,1,3)
plot(r.L3,r.max_kine,'o-')


ttt=1
