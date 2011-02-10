function vary_mass



%% Simulation parameters

% Precision of simulation (default = 10^-5, use 10^-7 for high precision)
p.maxError = 10^-7;

% Duration of simulation (s)
p.simDur    = 0.002;

% Time values to evaluate results(use 1000-5000)
p.t = linspace(0,p.simDur,1000);


visData = 0;
echoData = 0;
animate = 0;
indiv = 120;

p = get_params(indiv,p);

% Range of mass values
mass_min = p.dacMass * 0.7;
mass_max = p.dacMass * 1.3;

% Range of k_spring values
k_min = p.dacMass * 0.9;
k_max = p.dacMass * 1.1;

% Density of water (kg/m^3)
p.rho = 998;



%% Run for mass, normal L3

p = get_params(indiv,p);
p_start = p;


masses = linspace(mass_min,mass_max,10);

clrs = {'r-','g-','b-'};

mksize = 4;

for i = 1:length(masses)
    
    % Current L3 length
    p.dacMass = masses(i);
    p.dacI    = 0.28608 .* (( p.dacLen).^2) .* p.dacMass;
    
    %p.D = 2*p.D;
    
    %p.dacMass = p.dacMass*1.5;
    
    % Check geometry
    L = check_linkage(p,0);
    
    % Run the model with drag
    dD = run_sim(p,echoData);
    
    % Find index where max velocity is reached
    idx = find(dD.Dgamma == max(dD.Dgamma));
    
    % Store results
    rD.dacMass(i)    = p.dacMass;
    rD.gamma(i)      = max(dD.Dgamma);
    rD.E_drag(i)     = max(dD.E_drag);
    rD.theta_max(i)  = dD.gamma(idx)-dD.gamma(1);
    rD.t_peak(i)     = dD.t(dD.Dgamma==max(dD.Dgamma));
    rD.KT_min(i)     = min(dD.KT);
    rD.momentum(i)   = p.dacI .* max(dD.Dgamma);
    
    % Run the model without drag
    p.D = 0;
    dN = run_sim(p,echoData);
    
    % Find index where max velocity is reached
    idx = find(dD.Dgamma == max(dD.Dgamma));
    
    % Store results
    rN.dacMass(i)    = p.dacMass;
    rN.gamma(i)      = max(dN.Dgamma);
    rN.theta_max(i)  = dN.gamma(idx)-dN.gamma(1);
    rN.t_peak(i)     = dD.t(dN.Dgamma==max(dN.Dgamma));
    rN.KT_min(i)     = min(dN.KT);

    rN.al_pred(i)     = sqrt(p.kSpring./((p.dacI+p.waterI))) * ...
                             ((p.thetaRest-p.thetaStart));
    rN.momentum(i)   = p.dacI .* max(dN.Dgamma);
    %figure;plot(dN.t,p.thetaRest-dN.thetaIn)
    % Reset parameter values for next loop
    p = p_start;
    
    % Update status
    disp([num2str(i) ' of ' num2str(length(masses)) ' completed'])

end

figure;
% subplot(3,1,1)
% h = plot(1000.*rD.dacMass,1000.*rD.t_peak,'o-',...
%          1000.*rN.dacMass,1000.*rN.t_peak,'ro-');
% legend('w/drag','no drag')
% ylabel('time to peak speed (ms)')
% xlabel('Dac mass (g)')
% axis square
% 
% set(h(1),'MarkerSize',mksize)
% set(h(1),'MarkerFaceColor','b')
% set(h(2),'MarkerSize',mksize)
% set(h(2),'MarkerFaceColor','r')
% title('normal L3')

subplot(3,1,1)
h = plot(1000.*rD.dacMass,rD.theta_max.*180/pi,'o-',...
         1000.*rN.dacMass,rN.theta_max.*180/pi,'ro-');
legend('w/drag','no drag')
ylabel('Range of motion (deg)')
xlabel('Dac mass (g)')
axis square

set(h(1),'MarkerSize',mksize)
set(h(1),'MarkerFaceColor','b')
set(h(2),'MarkerSize',mksize)
set(h(2),'MarkerFaceColor','r')
title('normal L3')


subplot(3,1,2)
h = plot(1000.*rD.dacMass,rD.gamma.*(180/pi)./1000,'o-',...
         1000.*rN.dacMass,rN.gamma.*(180/pi)./1000,'ro-',...
         1000.*rN.dacMass,rN.al_pred.*(180/pi)./1000,'k-');
     
set(h(1),'MarkerSize',mksize)
set(h(1),'MarkerFaceColor','b')
set(h(2),'MarkerSize',mksize)
set(h(2),'MarkerFaceColor','r')

ylabel('max ang speed (deg/s)')
xlabel('Dac Mass (g)')
axis square

subplot(3,1,3)
h = plot(1000.*rD.dacMass,rD.momentum.*10^6,'o-',...
         1000.*rN.dacMass,rN.momentum.*10^6,'ro-');
     
set(h(1),'MarkerSize',mksize)
set(h(1),'MarkerFaceColor','b')
set(h(2),'MarkerSize',mksize)
set(h(2),'MarkerFaceColor','r')

ylabel('Angular momentum')
xlabel('Dac Mass (g)')
axis square

return

%% Run for mass, high L3

p = get_params(indiv,p);
p.L3 = p.L3 + 0.0003;
p_start = p;


masses = linspace(mass_min,mass_max,10);

clrs = {'r-','g-','b-'};

mksize = 3;

for i = 1:length(masses)
    
    % Current L3 length
    p.dacMass = masses(i);
    p.dacI    = 0.28608 .* (( p.dacLen).^2) .* p.dacMass;
    
    %p.D = 2*p.D;
    
    %p.dacMass = p.dacMass*1.5;
    
    % Check geometry
    L = check_linkage(p,0);
    
    % Run the model with drag
    dD = run_sim(p,echoData);
    
    % Find index where max velocity is reached
    idx = find(dD.Dgamma == max(dD.Dgamma));
    
    % Store results
    rD.dacMass(i)    = p.dacMass;
    rD.gamma(i)      = max(dD.Dgamma);
    rD.E_drag(i)     = max(dD.E_drag);
    rD.theta_max(i)  = dD.gamma(idx)-dD.gamma(1);
    rD.t_peak(i)     = dD.t(dD.Dgamma==max(dD.Dgamma));
    rD.KT_min(i)     = min(dD.KT);
    rD.momentum(i)   = p.dacI .* max(dD.Dgamma);
    
    % Run the model without drag
    p.D = 0;
    dN = run_sim(p,echoData);
    
    % Find index where max velocity is reached
    idx = find(dD.Dgamma == max(dD.Dgamma));
    
    % Store results
    rN.dacMass(i)    = p.dacMass;
    rN.gamma(i)      = max(dN.Dgamma);
    rN.theta_max(i)  = dN.gamma(idx)-dN.gamma(1);
    rN.t_peak(i)     = dD.t(dN.Dgamma==max(dN.Dgamma));
    rN.KT_min(i)     = min(dN.KT);

    rN.al_pred(i)     = sqrt(p.kSpring./((p.dacI+p.waterI))) * ...
                             ((p.thetaRest-p.thetaStart));
    rN.momentum(i)   = p.dacI .* max(dN.Dgamma);
    %figure;plot(dN.t,p.thetaRest-dN.thetaIn)
    % Reset parameter values for next loop
    p = p_start;
    
    % Update status
    disp([num2str(i) ' of ' num2str(length(masses)) ' completed'])

end


subplot(3,2,2)
h = plot(1000.*rD.dacMass,1000.*rD.t_peak,'o-',...
         1000.*rN.dacMass,1000.*rN.t_peak,'ro-');
legend('w/drag','no drag')
ylabel('time to peak speed (ms)')
xlabel('Dac mass (g)')
axis square

set(h(1),'MarkerSize',mksize)
set(h(1),'MarkerFaceColor','b')
set(h(2),'MarkerSize',mksize)
set(h(2),'MarkerFaceColor','r')
title('high L3')

subplot(3,2,4)
h = plot(1000.*rD.dacMass,rD.gamma.*(180/pi)./1000,'o-',...
         1000.*rN.dacMass,rN.gamma.*(180/pi)./1000,'ro-',...
         1000.*rN.dacMass,rN.al_pred.*(180/pi)./1000,'k-');
     
set(h(1),'MarkerSize',mksize)
set(h(1),'MarkerFaceColor','b')
set(h(2),'MarkerSize',mksize)
set(h(2),'MarkerFaceColor','r')

ylabel('max ang speed (deg/s)')
xlabel('Dac Mass (g)')
axis square

subplot(3,2,6)
h = plot(1000.*rD.dacMass,rD.momentum,'o-',...
         1000.*rN.dacMass,rN.momentum,'ro-');
     
set(h(1),'MarkerSize',mksize)
set(h(1),'MarkerFaceColor','b')
set(h(2),'MarkerSize',mksize)
set(h(2),'MarkerFaceColor','r')

ylabel('Angular momentum')
xlabel('Dac Mass (g)')
axis square


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
