function run_sims
% Performs calculations of relative flow from estimates of flow
% acceleration at the mouths of bluegill and bass.  Plots provide the basis
% for next-to-last figure in the Stewart & McHenry manuscript.



%% Define directories

dPath = '/Volumes/Docs/Projects/Relative velocity/lit_data';
zBaseM = '/Volumes/workgroup/Manuscripts/Relative flow/zMorphometrics_data';


%% Parameters

rho_water = 998;   % Density of water (kg m^-3)
time_for_flow = 0.015; % time at which to evaluate the relative flow (s) 


%% Load, spline fit, modify data from the literature

% Load structures blue and bass
load([dPath filesep 'lit_data.mat'])

% % Shift time of bass data
idx = find(bass.t>0.0004,1,'first');
bass.t = bass.t-bass.t(idx);

% Flip sign of U
bass.U = -(bass.U-bass.U(idx));
blue.U = -(blue.U-blue.U(1));


% Calculate velocities at the mouth
bass.Umouth = bass.U .* 4.59;
blue.Umouth = blue.U .* 3.56;

% Smoothing spline (for dynamic model)
blue.cs = csaps(blue.t,blue.Umouth);
bass.cs = csaps(bass.t,bass.Umouth);

% % Remove time before t=0
% idx = bass.t>0;
% bass.t = bass.t(idx);
% bass.U = bass.U(idx);
% bass.Umouth = bass.Umouth(idx);

clear idx 



% %% Data for impulse chamber
% 
% % Load stimulus data, stim
% load([dPath filesep 'stimulus_data.mat'])
% 
% % Define universal time vector
% t = linspace(-20e-3,0.9*min(stim.t(end,:)),size(stim.t,1))';
% 
% % Step through and interpolate data for t
% for i = 1:size(stim.t,2)
%     Ui(:,i) = interp1(stim.t(:,i),stim.U(:,i),t);
% end
% 
% % Store data
% impul.t = t;
% impul.U = mean(Ui,2)/1000;
% 
% clear t Ui i stim


%% Bluegill: Linear approximations for flow acceleration at the mouth

blue.idx = blue.t < 0.0183;
blue.c   = polyfit(blue.t(blue.idx),blue.Umouth(blue.idx),1);


%% Bass: Linear approximations for flow acceleration at the mouth

bass.idx = (bass.t < 0.021) & (bass.t > 0);
bass.c   = polyfit(bass.t(bass.idx),bass.Umouth(bass.idx),1);


%% Plot linear fits

symSize = 4;

figure;

% Bluegill: plot raw data
subplot(2,2,1:2)
h = plot(blue.t,blue.Umouth,'o');
set(h,'markerSize',symSize)
set(h,'markerFaceColor','b')
hold on

% Bluegill: linear fit
plot(blue.t(blue.idx),polyval(blue.c,blue.t(blue.idx)),'-')

% Bluegill: smoothing spline plot
fnplt(blue.cs,'-');


% Bass: plot raw data
h = plot(bass.t,bass.Umouth,'ro');
set(h,'markerSize',symSize)
set(h,'markerFaceColor','r')

% Bass: linear fit
plot(bass.t(bass.idx),polyval(bass.c,bass.t(bass.idx)),'r-')

% Bass: smoothign spline
fnplt(bass.cs,'r-');

% Title
title(['Bass accel = ' num2str(bass.c(1)) ...
       ' (m/s^2)  bluegill = ' num2str(blue.c(1)) ' (m/s^2)']) 

% Set limits, etc
xlim([0 0.06])
legend('blue',' ',' ','bass')
ylabel('Velocity (m/s)')
   
clear symSize h


%% Load specific gravity data

% Load pooled data, mP
load([zBaseM filesep 'body_metrics'])

% Identify those w/out an SB
idx = mP.sb_vol==0;

% Define mean body density (no SB)
rho_noSB = mP.rho_body(idx);

% Define mean body density (no SB)
rho_wSB = mP.rho_body(~idx);

% Calculate specific gravities
SG_noSB = rho_noSB / rho_water;
SG_wSB  = rho_wSB / rho_water;

clear mP rho_noSB rho_wSB


% %% Make predictions of relative flow
% 
% bass.Urel_noSB = time_for_flow .* bass.c(1) .* (1 - 1./SG_noSB);
% bass.Urel_wSB  = time_for_flow .* bass.c(1) .* (1 - 1./SG_wSB);
% blue.Urel_noSB = time_for_flow .* blue.c(1) .* (1 - 1./SG_noSB);
% blue.Urel_wSB  = time_for_flow .* blue.c(1) .* (1 - 1./SG_wSB);
% 
% X = [blue.Urel_noSB;       
%      bass.Urel_noSB; 
%      blue.Urel_wSB;
%      bass.Urel_wSB].*1000;
%  
% G = [ones(length(blue.Urel_noSB),1);     
%      2.*ones(length(bass.Urel_noSB),1);
%      3.*ones(length(blue.Urel_wSB),1);
%      4.*ones(length(bass.Urel_wSB),1)];
%  
% subplot(2,2,3)
% boxplot(X,G)
% ylabel('Relative velocity (mm/s)')
% xlabel('blue pre           bass pre           blue post          bass post')
% 
% disp('pre = pre inflation, post = post inflation')
% 


%% Make and plot predictions of relative flow

% Calculate relative flow 
bass.Urel_noSB = time_for_flow .* bass.c(1) .* (1 - 1./SG_noSB);
bass.Urel_wSB  = time_for_flow .* bass.c(1) .* (1 - 1./SG_wSB);
blue.Urel_noSB = time_for_flow .* blue.c(1) .* (1 - 1./SG_noSB);
blue.Urel_wSB  = time_for_flow .* blue.c(1) .* (1 - 1./SG_wSB);

% Organize dependent variable data for plot
X1 = [bass.Urel_noSB; 
      bass.Urel_wSB].*1000;
 
% Categorize data
G1 = [ones(length(bass.Urel_noSB),1);
     2.*ones(length(bass.Urel_wSB),1)];
 
% Organize dependent variable data for plot
X2 = [blue.Urel_noSB;       
     blue.Urel_wSB].*1000;
 
% Categorize data 
G2 = [ones(length(blue.Urel_noSB),1);     
     2.*ones(length(blue.Urel_wSB),1)];

% Plot for bass
subplot(2,2,3)
boxplot(X1,G1)
ylabel('Relative velocity (mm/s)')
xlabel('pre inflation                 post inflation')
title('Bass')
ylim([-120 0])

% Plot for bluegill
subplot(2,2,4)
boxplot(X2,G2)
ylabel('Relative velocity (mm/s)')
xlabel('pre inflation                   post inflation')
title('Bluegill')
ylim([-120 0])


%% Save data for dynamic simulations

save([dPath filesep 'lit_data_part2.mat'],'blue','bass')

return

%% Below is junk code that I'm not ready to throw away


% 
% %% Plot velocity and acceleration wrt position from the mouth
% 
% x = linspace(0,2.5,1000);
% 
% %4th order polynominal from Wainwright and Day (2007)...
% %that estimates fluid velocity in front
% %of predator's mouth as a function of
% %distance from mouth
% c1 = [0.098 -0.7 1.86 -2.19 1]; 
% 
% if 1
% figure;
% %subplot(2,1,1)
% plot(x,polyval(c1,x),'k')
% xlabel('position (gape diameters)')
%     ylabel('speed (relative to at the mouth)')
% 
% % subplot(2,1,2)
% % plot(x,polyval(polyder(c1),x),'r')
% 
% end
% 
% return
% 
% %% Dynamic simulations
% 
% 
% % Path to data
% zBaseM = '/Volumes/workgroup/Manuscripts/Relative flow/zMorphometrics_data';
% 
% % Density of water (kg m^-3)
% rho_water   = 998;     
% 
% % Load pooled data, mP
% load([zBaseM filesep 'body_metrics'])
% 
% % Identify those w/out a SB
% idx = mP.sb_vol==0;
% 
% % Define mean body density (no SB)
% rho_noSB = mean(mP.rho_body(idx));
% 
% % Define mean body density (no SB)
% rho_wSB = rho_water+0.*mean(mP.rho_body(~idx))
% 
% % Define I (with SB)
% I_wSB = mean(mP.I(~idx));
% 
% % Define body volume (no SB)
% V_noSB = mean(mP.Vbody(idx));
% 
% % Define body volume (with SB)
% V_wSB = mean(mP.Vbody(~idx));
% 
% % Define mass (with SB)
% M_wSB = mean(mP.Mbody(~idx));
% 
% % Mean COM with SB
% COM_wSB = mean(mP.COM(~idx,:));
% 
% % Distance fromCOM and tail tip
% tipLen = mean(mP.b_length(~idx))-mean(mP.COM(~idx,2));
% 
% % Mean position of level arm btn COV & COM with swim bladder
% L = [mean(mP.COV(~idx,1)-mP.COM(~idx,1)) ...
%     mean(mP.COV(~idx,2)-mP.COM(~idx,2)) ...
%     mean(mP.COV(~idx,1)-mP.COM(~idx,1))];
% 
% clear mP
% 
% 
% 
% % Bill's polynomial curve fit from Skorczewski et al (2010) that describes the
% % fluid speed at the predator's mouth during a strike as a function of time.
% % Time = 0 is begining of strike. 
% c1 = [-1.3211e10 3.4911e9 -3.9296e8 2.4207e7 -8.37e5 1.3109e4 8.2560 0.0087];
% 
% 
% %4th order polynominal from Wainwright and Day (2007)...
% %that estimates fluid velocity in front
% %of predator's mouth as a function of
% %distance from mouth
% c2 = [0.098 -0.7 1.86 -2.19 1]; 
% 
% % time vector
% t = linspace(0,.06,1000);
% 
% 
% % Simulation 1
% gapeDiameter = 2e-2;
% xInitial = 15e-3;
% thetaInitial = 45*(pi/180);
% 
% 
% % Initial conditions
% xPos = xInitial;
% U_body = 0;
% U_local = 0;
% U_relative = 0;
% theta = thetaInitial;
% alpha = 0;
% 
% % Calculate initial water acceleration
% U_mouth = - polyval(c1,t(2));
% U_local = U_mouth .* polyval(c2,xInitial./gapeDiameter);
% a_local = U_local./t(2);
% %t_a_water = t(2)/2;
% 
% clear U_mouth U_local
% 
% % Body coordinates
% bod    = [L(2).*sin(theta) L(2).*cos(theta) 0]; 
% a_bods = 0;
% a_locals = 0;
% 
% % Loop through time
% for i = 2:length(t)-1
%     
%     % Define change in time
%     dt = t(i)-t(i-1);
%     
%     % Pressure force and torque ___________________________________________
%     
%     Local pressure gradient
%     dpdx  = rho_water .* a_local;
%     
%     % Magnitude of pressure force
%     Fmag   = V_wSB.*dpdx;
%     
%     % Pressure force vector
%     F      = [0 Fmag 0]; 
% 
%     % Torque about COM
%     tau    = cross(F,bod); 
%        
%     
%     % Change in body position _____________________________________________
%     
%     % Translation state variables
%     a_bod       = dpdx/rho_wSB;
%     U_body(i,1) = U_body(i-1) + a_bod*dt;
%     xPos(i,1)   = xPos(i-1) + U_body(i).*dt;
%     a_bods(i,1) = a_bod;
%     
%     % Rotation state variables
%     omega       = tau(3)./I_wSB;
%     alpha(i,1)  = alpha(i-1) + omega*dt;
%     theta(i,1)  = theta(i-1) + alpha(i)*dt;
%     
%     % Update body coordinates with new theta
%     bod = [L(2)*sin(theta(i)) ...
%         L(2)*cos(theta(i)) ...
%         0];
%     
%     clear omega a_bod
%     
%     
%     % Calculate local and relative flow for current time __________________
%     
%     % Flow speed at the mouth (m/s)
%     U_mouth = - mean([ polyval(c1,t(i)) polyval(c1,t(i-1)) ]);
%     
%     % Flow speed at body position (m/s)
%     U_local(i) = U_mouth .* polyval(c2,mean([xPos(i-1)])./gapeDiameter);
%     
%     % Relative velocity at the COM
%     U_relative(i) = U_local(i) - U_body(i);
%     
%     clear U_mouth
%     
%     
%     % Calculate flow acceleration for next time ___________________________
%     
%     % Flow speed at the mouth (m/s)
%     U_mouth = -  mean([ polyval(c1,t(i)) polyval(c1,t(i+1)) ]);
%     
%     % Flow speed at body position (m/s)
%     U_next = U_mouth .* polyval(c2,xPos(i)./gapeDiameter);
%     
%     % Avg water acceleration near larva
%     a_local = (U_next - U_local(i))./(t(i+1)-t(i));
%  
%     a_locals(i) = a_local;
%     clear U_mouth U_next
%     
%    
%     if xPos <=0
%         return
%     end
%     
% %     subplot(2,1,1)
% % plot(1000.*t(1:i),1000.*(xPos(1:i)),'o-')
% % ylabel('Displacement (mm)')
% % subplot(2,1,2)
% % plot(1000.*t(1:i),U_body(1:i),'ro-',1000.*t(1:i),U_local(1:i),'o-')
% % xlabel('time (ms)')
% %pause
%     
% end
% 
% % Cut off last time value
% t = t(1:i);
% 
% figure;
% subplot(3,1,1)
% plot(1000.*t,- polyval(c1,t));
% ylabel('Velocity at the mouth')
% 
% subplot(3,1,2)
% plot(1000.*t,1000.*(xPos))
% ylabel('x-position (mm)')
% 
% subplot(3,1,3)
% plot(1000.*t,U_body,'r',1000.*t,U_local,'b',1000.*t,U_relative,'k-')
% ylabel('Velocity (m/s)')
% xlabel('time (ms)')
% legend('Body','Flow','Relative flow','Location','SouthWest')
% 
% 
% return
% % % Define simulation parameters
% % t           = [0:.0005:.08]';
% % flow_accel  = 1.6; % m s^-2
% % 
% % % Define pressure gradient
% % dpdx  = rho_water .* flow_accel;
% % 
% % % Simulation 1: translation, no SB
% % U_body1 = (1/rho_noSB) * dpdx .* t;
% % 
% % % Simulation 2: translation, with SB
% % U_body2 = (1/rho_wSB) * dpdx .* t;
% 
% % % Simulation 3: rotation_________________________________
% % theta  = pi/4;   % Body orientation wrt flow
% % Fmag   = V_wSB.*dpdx; % Magnitude of pressure force
% % bod    = [L(2).*sin(theta) L(2).*cos(theta) 0]; % Body coordinates
% % F      = [0 -Fmag 0]; % Force vector
% % tau    = cross(F,bod); % Torque about COM
% % omega  = tau(3)./I_wSB; % Rot Accelration at time = 0
% % alph   = 0; % Rot vel at time = 0
% % xAccel = F(2)/M_wSB; % Translations accel at time = 0
% % xVel   = 0; % x velocity at time = 0
% % xPos   = 0; % x position at time = 0
% 
% for i = 2:length(t)
%     % Time step
%     dt     = t(i)-t(i-1);
%     
%     % Current torque
%     tau    = cross(F,bod);
%     
%     % Rotation state variables
%     omega(i,1)  = tau(3)./I_wSB;
%     alph(i,1)   = alph(end) + omega(end)*dt;
%     theta(i,1)  = theta(end) + alph(end)*dt;
%     
%     % Translation state variables
%     xAccel(i,1) = F(2)/M_wSB;
%     xVel(i,1)   = xVel(end) + xAccel(end).*dt;
%     xPos(i,1)   = xPos(end) + xVel(end).*dt;
%     
%     % Update coordinates
%     bod = [L(2)*sin(theta(i)) ...
%         L(2)*cos(theta(i)) ...
%         0];
% end
% 
% % Plot kinematics
% figure;
% subplot(4,1,1)
% plot(t.*1000,flow_accel.*t.*1000,'b',...
%     t.*1000,U_body1.*1000,'r',...
%     t.*1000,U_body2.*1000,'r--')
% xlabel('t (ms)')
% ylabel('U (mm/s)')
% 
% subplot(4,1,2)
% plot(t.*1000,theta.*(180/pi))
% xlabel('t (ms)')
% ylabel('theta (deg)')
% 
% subplot(4,1,3)
% plot(t.*1000,xPos.*1000)
% xlabel('t (ms)')
% ylabel('Position (mm)')
% 
% subplot(4,1,4)
% plot(t.*1000,flow_accel.*t.*1000-U_body2.*1000,'r');
% hold on
% plot(t.*1000,alph.*tipLen.*1000,'b')
% hold off
% ylabel('velocity (mm/s)')
% xlabel('t (ms)')
% legend('Relative flow','rot of tail tip','Location','NorthWest')
% 
% clear idx rho_noSB rho_wSB I_wSB V_noSB V_wSB COM_wSB L t flow_accel
% clear dpdx U_body1 Ubody2 bodyS i S pForce omega alpha Fmag bod F tau
% clear xAccel xVel xPos dt