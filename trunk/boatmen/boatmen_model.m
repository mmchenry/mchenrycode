function r = boatmen_model(pIn,Cd_app,Cd_body)
% ODE model of boatmen hydrodynamics.
% This is a modified version of Blake, 1985
% pIn - Inout parameters required to run a simulation
% r   - structure of simulaton results
% sim_mode - ('optimize','default')


%% Input parameter values

% Body length (m) 
L = pIn.body_len;

% Initial body speed (m/s)
U0 = pIn.initial_speed;

% Angle of appendage at start of power stroke (rad)
%gamma1  = pIn.ang_start;

% Beat period (s)
%P = pIn.beat_period;

% Length of the paddle (m)
pLen = pIn.app_len;

% Amplitude of appendage speed (m/s)
spd_A = pIn.spd_A;

% Initial speed of power stroke (m/s)
spd_0 = pIn.spd_0;

% Phase of speed function (s)
spd_phs   = pIn.spd_phs;

% Period of speed function (s)
spd_P     = pIn.spd_P;

% Period of angle function (s)
ang_P     = pIn.ang_P;

% Amplitude of angle (rad)
ang_amp   = pIn.ang_amp;

% Start angle (rad)
ang_start = pIn.ang_start;

% Body mass (kg)
%mass_body = 23e-6;
mass_body = pIn.body_mass;

% Period of the return stroke
rtrn_P = pIn.rtrn_P;

% Number of strokes
num_strokes = pIn.num_strokes;

% Parameter defaults
if nargin < 3  
    % Drag coefficient of the body
    Cd_body = 1.07;
    
    if nargin < 2
        % Drag coefficient of the appendage
        Cd_app = 1.1;
         
    end
end

clear pIn


%% General & calculated parameter values

% Body density (kg m^-3)
%rho_body  = 1030;

% Water density (kg m^-3)
rho_water = 1000;

% Viscosity of water (Pa s)
mu = 0.001;

% Projected area of body (m^2) -- calculation assumes same shape as Blake
S = pi*(0.2*L)^2;

% Time to stop evaluation (s)
%t_stop = 100e-3;

% Initial body position (m)
X0 = 0;



% Height of appendage (m) -- calculation assumes same shape as Blake
h = 0.23.*L;


%% Scale parameters to non-dimensional units

% Define scaling constants
sL = L;
sT = 10^-3;
sM = mass_body*10^6;
sF = sM .* sL ./ sT^2;

% Declare global variables
%global p

% Scale all parameters
p.rho_water = rho_water   /sM *sL^3;
%p.rho_body  = rho_body    /sM *sL^3;
p.mass_body = mass_body   /sM;
cp.mu        = mu          /(sF/sL^2) /sT;
p.L         = L           /sL;
%p.tspan     = [0 t_stop]  ./sT;
p.U0        = U0          /sL *sT;
p.X0        = X0          /sL;
%p.gamma1    = gamma1;
%p.P         = P           /sT;
p.Cd_app    = Cd_app;
p.Cd_body   = Cd_body;
p.h         = h           /sL;
p.pLen      = pLen        /sL;
p.S         = S           /(sL^2);
p.spd_A     = spd_A       /sL *sT;
p.spd_phs   = spd_phs     /sT;
p.spd_P     = spd_P       /sT;
p.ang_P     = ang_P       /sT;
p.spd_0     = spd_0       /sL *sT;
p.ang_amp   = ang_amp;
p.ang_start = ang_start;
p.rtrn_P    = rtrn_P      /sT; 
p.num_strokes = num_strokes;

%clear mu L k t_stop rho_body rho_water S Cd_app Cd_body


%% Calculated & solver parameters

% Body volume (m^3)
%p.vol_body = (4/3) * pi * (p.L/2)^3;

% Solver options
options    = odeset('RelTol',1e-7);


%% Run numerical simulation 

[t,X] = solver(p,options);
%[t,X] = ode45(@gov_eqn,p.tspan,[p.X0; p.U0]);


%% Store results (convert values back to SI units)


% Appendage angle and speed 
gama  = ang_func(t,p.ang_amp,p.spd_P,p.ang_start,p.spd_phs,p.rtrn_P);
v_n   = spd_func(t,p.spd_A,p.spd_P,p.spd_phs,p.rtrn_P,p.spd_0);

% Thrust
T = app_thrust(gama,v_n,p,X(:,2));
D = body_drag(X(:,2),p);

%if strcmp(sim_mode,'default')
    
    % return a full set of result variables
    r.v_n    = v_n     .*sL ./sT;
    r.gamma  = gama;
    r.thrust = T       .*sF;
    r.drag   = D       .*sF;
    r.t      = t       .*sT;
    r.x      = X(:,1)  .*sL;
    r.U      = X(:,2)  .*sL  ./sT;
    
%     % Report results    
%     ReMax_body   = L .* max(abs(r.U)) .* rho_water ./ mu;
%     ReMax_app    = max(abs(r.v_n)) .* pLen .* rho_water ./ mu;
%     
%     disp(' ')
%     disp([' Maximum Re of the body = ' num2str(ReMax_body)]);
%     disp(' ')
%     disp([' Maximum Re of the appendage = ' num2str(ReMax_app)]);
%     disp(' ')
%     
    % Visualize results
    if 0
        % Plot
        figure;
        subplot(3,1,1)
        [ax,h1,h2] = plotyy(r.t.*1000,r.gamma,r.t.*1000,r.v_n.*1000);
        ylabel(ax(1),'gamma (rad)')
        ylabel(ax(2),'v_n (mm/s)')
        grid on
        
        subplot(3,1,2)
        plot(r.t.*1000,r.thrust.*10^6,'g',r.t.*1000,r.drag.*10^6,'-')
        legend('thrust','drag')
        ylabel('Force (micro N)')
        grid on
        
        subplot(3,1,3)
        plot(r.t.*1000,r.U.*1000,'-')
        grid on
        ylabel('U')
        
    end
    
% elseif strcmp(sim_mode,'optimize')
%     
%     % Calculate speed and time in SI units
%     U = X(:,2)  .*sL  ./sT;
%     tU = t       .*sT;
%     
%     % Calculate mean speed after initial stroke
%     idx = tU > (spd_P/2 + rtrn_P);
%     
%     r = mean(U(idx));
%     
% else
%     
%     error('Requested sim_mode not an option')
%     
% end

clear U t X T D




% subplot(3,1,3)
% [ax2,h3,h4] = plotyy(r.t.*1000,r.U.*1000,r.t.*1000,r.alpha.*l_tip.*1000);
% %legend('body','appendage')
% xlabel('time (ms)')
% ylabel(ax2(1),'Body speed (mm/s)')
% ylabel(ax2(2),'Appendage speed (mm/s)')
% ylim(ax2(1),[0 800])
% ylim(ax2(2),[-800 0])
% grid on

end

function [t,X] = solver(p,options)

tspan = [0 p.num_strokes*(p.spd_P+p.rtrn_P)];

[t,X] = ode45(@gov_eqn,tspan,[p.X0; p.U0],options);

    function dX = gov_eqn(t,X)
        % ODE of the dynamics of the system
        
        % Body position
        pos = X(1);
        
        % Body speed
        U   = X(2);
        
        alph = 0;
        
        gama  = ang_func(t,p.ang_amp,p.spd_P,p.ang_start,p.spd_phs,p.rtrn_P);
        v_n   = spd_func(t,p.spd_A,p.spd_P,p.spd_phs,p.rtrn_P,p.spd_0);
        
        % Appendage position (gama) and speed (alph)
        %[gama,alph] = app_kine(t,p.gamma1,p.gamma2,p.P);
        
        % Thrust
        thrust = app_thrust(gama,v_n,p,U);

        % Drag
        drag = body_drag(U,p);
        
        % Define output: speed
        dX(1,1) = U;
        
        % Define output: acceleration
        dX(2,1) = (thrust + drag) ./ p.mass_body;
        
    end

end

function thrust = app_thrust(gama,v_n,p,U)
% Calculate thrust generated by appendage

v_app = - v_n + U.*cos(gama-pi/2);

thrust = - cos(gama-pi/2)*....
                    0.5*p.rho_water.*p.Cd_app.*...
                    p.h.*p.pLen.*v_app.*abs(v_n);
     %TODO: subtract velocity from forward movement
     
% Adjust direction of thrust
%thrust = (gama-pi/2)./abs((gama-pi/2)) .* thrust;

end
  
function D = body_drag(U,p) 
% Calculates drag force on the body

D = -0.5 .* p.S .* p.rho_water .* p.Cd_body .* abs(U) .* U;

end

function y = spd_func(t,A,P,phs,rtrn_P,spd_0)
% Function that defines the speed of the power stroke

%phs = 0;
% Make sure t is a column vector
if size(t,2) > size(t,1)
    t = t';
end

% Define which stroke number we are within
stroke_num = floor(t/(P+rtrn_P))+1;

% Time adjusted to stroke cycle
t_ad = (t - ((stroke_num-1) * (P+rtrn_P)));

% Index of power stroke times
idx = t_ad <  (P - 0*phs);

% Power stroke velocity values
y(idx,1) = spd_0 + A.*sin(pi.*(t_ad(idx)-phs)./(1.5*P)).^2;

% Zero values for recovery stroke
y(~idx,1) = zeros(size(sum(~idx),1),size(sum(~idx),2));

end

function y = ang_func(t,A,P,ang_start,phs,rtrn_P)
% Defines angle of power stroke
%phs=0;
% Make sure t is a column vector
if size(t,2) > size(t,1)
    t = t';
end

% Define which stroke number we are within
stroke_num = floor(t/(P+rtrn_P))+1;

% Time adjusted to stroke cycle
t_ad = (t - ((stroke_num-1) * (P+rtrn_P)));

% Index of power stroke times
idx = t_ad <  (P - 0*phs);

% Angle for power stroke
y(idx,1) = ang_start + A.*sin(pi.*(t_ad(idx))./(2.*(1.1*P))).^2;

% Zero values for recovery stroke
y(~idx,1) = zeros(size(sum(~idx),1),size(sum(~idx),2));

end
