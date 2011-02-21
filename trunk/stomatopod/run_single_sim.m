function run_single_sim
% Runs the model and visualizes the resutls for a single individual
% Similar to model_test, but doesn't rerun for all indivduals

% Indivudal to simulate
indiv = 120;



%% Simulation parameters

% Precision of simulation (default = 10^-5, use 10^-7 for high precision)
p.maxError = 10^-7;

% Duration of simulation (s)
p.simDur    = 0.02;

% Time values to evaluate results(use 1000-5000)
p.t = linspace(0,p.simDur,1000);


%% General parameter values
% Fixed for all individuals

% water I (kg m^2)
%p.waterI      = 1e-8;

% 3rd moment of area of dactyl (m^5) (for 12 mm long dactyl, r = 3 mm)
%p.dacA         = 0.25*(3e-3)*(12e-3)^4;

% Density of water (kg/m^3)
p.rho = 998;

%TODO: Adjust these to be custom for each individual


%% Run simulation


% Individual-specific parameter values
p = get_params(indiv,p);

% Knock out drag
%p.D = 0;

%p.L3 = p.L3 + 0.0005;
%p.L3 = p.L3 - .0001;
%p.waterI = 0;

%dL = 1e-4;
%p.dacMass  = p.dacMass * 0.4;
%p.dacI    = 0.28608 .* (( p.dacLen).^2) .* p.dacMass;
%p.D = 3*p.D;
p.rel_tol = 1e-7;
% Run simulation
[d,result] = run_sim(p,0);

L = check_linkage(p,0,10^6);
minKT = min(L.KT);




%% Visualize simulation results (v. 1)
figure;

subplot(4,1,1)
plot(d.t.*1000,d.theta.*180/pi,'r')
ylabel('theta (deg)')
xlabel('time (ms)')
grid on

subplot(4,1,2)
plot(d.t.*1000,(d.gamma)*180/pi,'b')
ylabel('gamma (deg)')
xlabel('time (ms)')
grid on

% 
% subplot(3,1,2)
% plot(d.t.*1000,(d.dacAngSpd.*180/pi)./1000)
% ylabel('SB angular spd (deg/ms)')
% xlabel('time (ms)')
% grid on

subplot(4,1,3:4)
E_tot = d.E_kin + d.E_drag + d.E_elastic;
%E_tot = d.E_kin + d.E_elastic;
plot(d.t.*1000,d.E_elastic.*1000,'r',d.t.*1000,d.E_drag.*1000,'b',...
    d.t.*1000,d.E_kin.*1000,'g',d.t.*1000,E_tot.*1000,'k--')
%plot(d.t.*1000,d.E_elastic.*1000,'r',d.t.*1000,d.E_kin.*1000,'g',d.t.*1000,E_tot.*1000,'k--')
ylabel('Energy (mJ)')
xlabel('time (ms)')
%grid on
legend('elastic','drag','kinetic','total','Location','West');
title(['min KT = ' num2str(minKT)])
%xlim([0 2])

% % Power graph
% subplot(2,1,2)
% plot(d.t(2:end).*1000,diff(d.E_elastic)./diff(d.t),'r',...
%     d.t(2:end).*1000,diff(d.E_drag)./diff(d.t),'b',...
%     d.t(2:end).*1000,diff(d.E_kin)./diff(d.t),'g')
% ylabel('Power (J/s)')
% xlabel('time (ms)')
% grid on


clear E_tot


return

%% Visualize simulation results (v. 2)
figure;

subplot(4,1,1)
plot(d.t.*1000,(d.Dgamma)*180/pi,'b')
ylabel('Dgamma (deg/s)')
xlabel('time (ms)')
grid on
xlim([0 2])
title(num2str(minKT))

subplot(4,1,2)
plot(d.gamma,(d.Dgamma)*180/pi,'b')
xlabel('gamma (deg)')
ylabel('Dgamma^2 (deg/s)')
grid on
%ylim([0 15e7])

% 
% subplot(3,1,2)
% plot(d.t.*1000,(d.dacAngSpd.*180/pi)./1000)
% ylabel('SB angular spd (deg/ms)')
% xlabel('time (ms)')
% grid on

subplot(4,1,3)
E_tot = d.E_kin + d.E_drag + d.E_elastic;
%E_tot = d.E_kin + d.E_elastic;
plot(d.t.*1000,d.E_elastic.*1000,'r',d.t.*1000,d.E_drag.*1000,'b',...
    d.t.*1000,d.E_kin.*1000,'g',d.t.*1000,E_tot.*1000,'k--')
%plot(d.t.*1000,d.E_elastic.*1000,'r',d.t.*1000,d.E_kin.*1000,'g',d.t.*1000,E_tot.*1000,'k--')
ylabel('Energy (mJ)')
xlabel('time (ms)')
%grid on
legend('elastic','drag','kinetic','total','Location','West');
title(['min KT = ' num2str(minKT)])
xlim([0 2])
%xlim([0 2])

% Power graph
subplot(4,1,4)
plot(d.t(2:end).*1000,diff(d.E_elastic)./diff(d.t),'r',...
    d.t(2:end).*1000,diff(d.E_drag)./diff(d.t),'b',...
    d.t(2:end).*1000,diff(d.E_kin)./diff(d.t),'g')
ylabel('Power (J/s)')
xlabel('time (ms)')
grid on
xlim([0 2])

return


%% Animate results

f = figure;
set(f,'DoubleBuffer','on');
pause(.1)
for i = 1:size(d.A,1)
    h = plot([d.A(i,1) d.B(i,1)],[d.A(i,2) d.B(i,2)],'b-',...
        [d.C(i,1) d.B(i,1)],[d.C(i,2) d.B(i,2)],'b-',...
        [d.C(i,1) d.D(i,1)],[d.C(i,2) d.D(i,2)],'b-');
    xlim([0 6e-3])
    axis equal
    pause
end


%TODO: calc amd compare momentum values
%TODO: Set up to run all
%TODO: 




