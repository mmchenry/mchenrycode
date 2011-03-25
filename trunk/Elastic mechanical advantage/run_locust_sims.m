function run_locust_sims

% Marker size (in plots)
mSize = 3;


%% Vary in-lever length

% Default parameters
p = default_param;

% Set values for in-levers
%in_levers = p.hIn .* 10.^[-0.5:.3:0.5];
in_levers = [.05:.01:0.1].*0.0173;

% Remove fluid forces and gravity
p.g   = 0;
p.rho = 0;
    
%p.L_s = 1.2e-3;
p.k = 1e4;

for i = 1:length(in_levers)
    
    % Set current in-lever
    p.hIn = in_levers(i);
    
    [d,p] = locust_model(p,0,0);
    
    max_v(i) = max(d.v);
    max_t(i) = d.t(end-1);
    
    %EMA(i) = max(d.EMA);
    EMA(i) = p.hIn/p.hOut;
    hit_max(i) = d.hit_max;
    
    pre_max_t(i) = pi/EMA(i) * sqrt(p.m/p.k);
    %pre_max_t(i) = pi/(p.hIn/p.hOut) * sqrt(p.m/p.k);
    
    pre_max_v(i) = d.L_range * sqrt(p.k/p.m);
    clear d
    
end



f1 = figure;

% Max velocity plot ---------------------------------
subplot(1,2,1)
plot(EMA(hit_max),max_v(hit_max),'xr')
hold on
h = plot(EMA(~hit_max),max_v(~hit_max),'ro');
set(h,'markerEdgeColor','none')
set(h,'markerFaceColor','r')
set(h,'markerSize',mSize)
hold on
plot(EMA(~hit_max),pre_max_v(~hit_max),'k-')
xlabel('MA')
ylabel('max vel')
%legend('partial stretch','full stretch','Location','Southeast')
axis square

% Time to max plot ---------------------------------
c = polyfit(log10(EMA),log10(max_t),1);

subplot(1,2,2)
h = plot(log10(EMA),log10(max_t),'ro');
set(h,'markerEdgeColor','none')
set(h,'markerFaceColor','r')
set(h,'markerSize',mSize)
hold on
plot(log10(EMA),c(1).*log10(EMA)+c(2),'r-')
axis square

xlabel('log10(MA)')
ylabel('log10(lambda)')

clear in_levers pre_max_v pre_max_t EMA p hit_max max_v max_t


%% Vary out-lever length

% Default parameters
p = default_param;

% Set values for in-levers
%out_levers = p.hIn .* 10.^[-0.5:.3:0.5];
MAs = [.01:.01:0.06]+.005;

% Remove fluid forces and gravity
p.g   = 0;
p.rho = 0;
    
%p.L_s = 1.2e-3;
p.k = 1e4;

for i = 1:length(MAs)
    
    % Set current in-lever
    p.leg_len = sqrt(3*(p.hIn/MAs(i))^2);
    
    [d,p] = locust_model(p,0,0);
    
    p.hOut*1000
    max_v(i) = max(d.v);
    max_t(i) = d.t(end-1);
    
    %EMA(i) = max(d.EMA);
    EMA(i) = p.hIn/p.hOut;
    hit_max(i) = d.hit_max;
    
    pre_max_t(i) = pi/EMA(i) * sqrt(p.m/p.k);
    %pre_max_t(i) = pi/(p.hIn/p.hOut) * sqrt(p.m/p.k);
    
    pre_max_v(i) = d.L_range * sqrt(p.k/p.m);
    clear d
    
end

p       = default_param;
[d,p]   = locust_model(p,0,0);
MA_real = p.hIn/p.hOut;


figure(f1);
%figure;

subplot(1,2,1)
plot(EMA(hit_max),max_v(hit_max),'x')
hold on
h = plot(EMA(~hit_max),max_v(~hit_max),'bo');
set(h,'markerEdgeColor','none')
set(h,'markerFaceColor','b')
set(h,'markerSize',mSize)
set(gca,'YTick',[31:1:34])

ylim([30 35])
xlim([0 .11])
plot(EMA(~hit_max),pre_max_v(~hit_max),'k-')
plot([MA_real MA_real],ylim,'k--')
hold off
xlabel('MA')
ylabel('max vel')
%legend('partial stretch','full stretch','Location','Southeast')
axis square


% Time to max plot ---------------------------------
c = polyfit(log10(EMA),log10(max_t),1);

subplot(1,2,2)
h = plot(log10(EMA),log10(max_t),'o');
set(h,'markerEdgeColor','none')
set(h,'markerFaceColor','b')
set(h,'markerSize',mSize)

hold on
plot(log10(EMA),c(1).*log10(EMA)+c(2),'-')
axis square

xlim([-2 -0.8])
xlabel('log10(MA)')
ylabel('log10(lambda)')



return

% Scaling of period
c = polyfit(log10(out_levers(~hit_max)./p.hOut),log10(max_t(~hit_max)),1);


f2 = figure;
subplot(1,2,1)
plot(log10(out_levers(~hit_max)./p.hOut),log10(max_t(~hit_max)),'b-o',...
     log10(out_levers(~hit_max)./p.hOut),...
     log10(c(2).*(out_levers(~hit_max)./p.hOut).^c(1)),'g--')
title(['fit = ' num2str(c(1)) ])




% subplot(1,2,2)
% plot(EMA(hit_max),1000.*max_t(hit_max),'o',...
%      EMA(~hit_max),1000.*max_t(~hit_max),'ro')
% hold on
% plot([MA_real MA_real],ylim,'k--')
% plot(EMA(~hit_max),1000.*pre_max_t(~hit_max),'k-')
% hold off
% xlabel('MA')
% ylabel('time to max vel (ms)')
% axis square

% figure
% subplot(1,2,1)
% plot(pre_max_t(~hit_max),max_t(~hit_max),'o',...
%      pre_max_t(~hit_max),pre_max_t(~hit_max),'-')
% axis equal
% xlabel('Analytical')
% ylabel('Numerical')
% title('Period')

% subplot(1,2,2)
% plot(pre_max_v(~hit_max),max_v(~hit_max),'o',...
%      pre_max_v(~hit_max),pre_max_v(~hit_max),'-')
% axis equal
% xlabel('Analytical')
% ylabel('Numerical')
% title('V_m_a_x')

clear in_levers pre_max_v pre_max_t EMA p out_levers



return


%% Simulated experimental manipulation

p = default_param;

disp('Cutting the leg ------------------------------')
disp(' ')

p.m = 10e-6;
p.leg_len = (10/21)*p.leg_len;
p.L_s = 0.97e-3;
p.k = .73e4;

d = locust_model(p,0,0);

disp(' ')
disp('Adding a weight ------------------------------')

p.m = 210e-6;

d = locust_model(p,0,0);



function p = calc_params(p)
% Calculated parameters

p.xJoint = p.L_s-p.hIn.*cos(p.theta0);

% Moment of inertia (kg m^2)
%p.I = (p.m/12) * (p.leg_len^2 + p.leg_dia^2);
p.I = (p.m/3) * (p.leg_len^2);

% Moment of water inertia (kg m^2)
p.I_water = pi/3 * p.rho * (p.leg_dia/2)^2 * p.leg_len^3;

% Length of out-lever (m)
p.hOut = sqrt(p.I/p.m);

% Drag product
p.drag_prod = (pi/4) * p.Cd * p.rho * (p.leg_dia/2) * p.leg_len^4;



