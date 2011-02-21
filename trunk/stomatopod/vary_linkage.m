function vary_linkage



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
indiv = 120;

p = get_params(indiv,p);

L3_min = p.L3 - 1e-4;
L3_max = p.L3 + 10e-4;



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


L3s = linspace(L3_min,L3_max,10);
L3s = L3s(end:-1:1);

clrs = {'r-','g-','b-'};

mksize = 4;

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
    idx = find(dD.Dgamma == max(dD.Dgamma));
    
    % Store results
    rD.L3(i)         = p.L3;
    rD.gamma(i)      = max(dD.Dgamma);
    rD.E_drag(i)     = max(dD.E_drag);
    rD.theta_max(i)  = dD.gamma(idx)-dD.gamma(1);
    rD.t_peak(i)     = dD.t(dD.Dgamma==max(dD.Dgamma));
    rD.KT_min(i)     = min(dD.KT);
    rD.efficiency(i) = 100*max(dD.E_kin)./max(dD.E_elastic);
    
    % Run the model without drag
    p.D = 0;
    dN = run_sim(p,echoData);
    
    % Find index where max velocity is reached
    idx = find(dD.Dgamma == max(dD.Dgamma));
    
    % Store results
    rN.L3(i)         = p.L3;
    rN.gamma(i)      = max(dN.Dgamma);
    
    rN.theta_max(i)  = dN.gamma(idx)-dN.gamma(1);
    rN.t_peak(i)     = dD.t(dN.Dgamma==max(dN.Dgamma));
    rN.KT_min(i)     = min(dN.KT);

    rN.al_pred(i)     = sqrt(p.kSpring./((p.dacI+p.waterI))) * ...
                             ((p.thetaRest-p.thetaStart));
    rN.efficiency(i) = 100*max(dN.E_kin)./max(dN.E_elastic);
    %figure;plot(dN.t,p.thetaRest-dN.thetaIn)
    % Reset parameter values for next loop
    p = p_start;
    
    % Update status
    disp([num2str(i) ' of ' num2str(length(L3s)) ' completed'])

end


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
