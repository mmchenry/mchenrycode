function create_data

% Test our implementation of Lighthill's kinematic analysis by generating
% head kinematics


%% Parameter values

% Preliminary
headLength  = 2 ;       % Head length (cm)
headWidth   = .75;      % Head width (=2b from Lighthill, cm)
f           = 4;        % Beat freq (Hz)
framerate   = 250;      % Frame rate (fps)
numCycles   = 40;       % Number of cycles considered

% Position
% lat_amp     = headLength/10 *10^-2;    % Lateral amplitude (m)
% for_amp     = headLength/50 *10^-2;    % Fore-aft amplitude (m)
alpha_amp   = 10.*(pi/180);              % Yaw amplitude (rad)

% Velocity
U_inf       = 4.5;                % Flow tank speed (cm/s)
V_amp       = 5;                  % Amp of lateral velocity (cm/s)
U_amp       = 1;                % Amp of fore-aft velocity (cm/s)
omega_amp   = alpha_amp*(2*pi*f); % Amp of yaw rate (rad/s)
phi         = pi/2;             % Phase shift of yaw (rad)

% Output
visData         = 1;    % Provide timeseries graphs (1 or 0) 
makeMovie       = 0;    % Create animation of the kinematics (1 or 0)
interactiveMode = 0;    % Allows user to step trough animation (1 or 0)
noise_amp       = 0.001;% Proportion of noise added to coordinates

% Pressure
K   = 3;    % Numerical constant (Lighthill used K=3)
k   = 1.6;  % Shape factor (Lighthill used k=1.6)
rho = 998;  % Water density (kg/m^3)

disp(' ')
disp(['k ratio = ' num2str(U_inf*alpha_amp/V_amp)])
disp(' ')

%% Kinematic equations

t       = linspace(0,numCycles.*(1/f),numCycles.*(1/f)*framerate)';

% Position data.  Uses a different equation for alpha so that phi is
% correct wrt Lighthill (phase in alpha wrt V);
alpha_f = alpha_amp.*cos(2*pi*f.*t + phi); 

% Velocity data
V       = V_amp      .*cos(2*pi*f.*t);
U       = U_amp      .*cos((2*pi*f/2).*t);
omega   = -omega_amp .*sin(2*pi*f.*t + phi);

% Acceleration data
A       = -V_amp.*(2*pi*f)^2 .* sin(2*pi*f.*t);


%% Pressure calculation
dP  = headWidth*(1+k)*rho.*A - headWidth*K*rho.*(U+U_inf).*omega;


%% Display timeseries of data

if visData
   figure;
   
   subplot(6,1,1)
        plot(t,U,'.')
        grid on
        ylabel('U');
   subplot(6,1,2)
        plot(t,V,'.')
        grid on
        ylabel('V');
   subplot(6,1,3)
        plot(t,A,'.')
        grid on
        ylabel('A');
   subplot(6,1,4)
        plot(t,alpha_f,'.')
        grid on
        ylabel('alpha');
   subplot(6,1,5)
        plot(t,omega.*(180/pi),'.')
        grid on
        ylabel('omega');
   subplot(6,1,6)
        plot(t,dP,'.')
        grid on
        ylabel('dP');
        xlabel('time (s)')
end


%% Kinematic simulation

dt  = mean(diff(t));

alpha(1) = alpha_f(1);
x0(1)    = 0;
y0(1)    = 0;
x1(1)    = x0(1) - headLength.*cos(alpha(1));
y1(1)    = y0(1) + headLength.*sin(alpha(1));

if visData 
    figure
    set(gcf,'DoubleBuffer','on')
end

for i = 2:length(t)
    
    % Positon change of P0 in local system
    dX0    = mean([U(i-1) U(i)])*dt;
    dY0    = mean([V(i-1) V(i)])*dt; 
    
    % Positon changes of P0 in inertial system 
    dy0    = dY0*cos(alpha(i-1));
    dx0    = dX0*cos(alpha(i-1));
    
    % Change in alpha
    dAlpha = mean([omega(i-1) omega(i)]).* dt;
    
    % Store coordinates 
    x0(i)    = x0(i-1) + dx0;
    y0(i)    = y0(i-1) + dy0;
    alpha(i) = alpha(i-1)+dAlpha;
    x1(i)    = x0(i) - headLength.*cos(alpha(i));
    y1(i)    = y0(i) + headLength.*sin(alpha(i));
        
    clear d0 d1 dX1 dY1 dX0 dY0 orgn S dAlpha
end

% Clean up
clear dt
if visData 
    close
end


%% Add noise

x0  = x0 + range(x0).*noise_amp.*rand(size(x0));
y0  = y0 + range(y0).*noise_amp.*rand(size(y0));
x1  = x1 + range(x1).*noise_amp.*rand(size(x1));
y1  = y1 + range(y1).*noise_amp.*rand(size(y1));


%% Visualize kinematics

axRange = [-1.25*headLength 0.25*headLength -0.75*headLength 0.75*headLength];

if makeMovie
    % Animation
    fig = figure;
    set(fig,'DoubleBuffer','on');
    subplot(2,1,1);subplot(2,1,2);axis square

    if interactiveMode
        j = 1;
        tInterval = (1/f);   %Interval over which to view timeseries
        disp(' ')
        disp('right arrow - advance frame')
        disp('left arrow - back one frame')
        disp('return - quits');
        disp(' ')
    else
        tInterval = 2*(1/f); %Interval over which to view timeseries
    end

    for i=1:length(x0)
        
        if ~interactiveMode
            j = i;
        end
        
        subplot(2,1,1)
            plot(t,alpha./alpha_amp,'r-',t,V./V_amp,'b--',...
                [t(j) t(j)],[-1 1],'k');
            
            title('Normalized alpha (red)   V (blue)')
            set(gca,'XLim',[t(j)-tInterval/2 t(j)+tInterval/2]);
            grid on
            %legend('alpha','V','Location','NorthWest');
        subplot(2,1,2)
            plot([x0(j) x1(j)],[y0(j) y1(j)],'-',x0(j),y0(j),'o')
            hold on
            h = plot(x0(1:j),y0(1:j),'-');
            set(h,'Color',.5'*[1 1 1]);
            hold off
        
        %axis square
        axis(axRange);
        c = get(fig,'Children');
        set(c(1),'Position',[.2 .1 .6 .6])
        set(c(2),'Position',[.13 .75 .775 .15])
        grid on
        
        if interactiveMode
            
            [tX,tY,but] = ginput(1);
            
            if isempty(but) % Return
                return
            elseif but==29  % Right arrow
                j = min([j+1 length(x0)]);
            elseif but==28  % Left arrow
                j = max([j-1 1]); 
            end
            
        else
            pause(.3)
        end
        
        clear c
    end
end

if visData
    %Plots
    figure;
    
    plot(x0,y0);
    hold on
    plot(x1,y1,'r-')

    axis(axRange);
    grid on
    legend('point 0','point 1')
end

%disp('

%% Save data

d.xCtr = x0';
d.yCtr = y0';
d.xTip = x1';
d.yTip = y1';
d.framenum = [1:length(x0)]';

p.calconst = 1;
p.framerate = framerate;
p.units = 'cm';
p.motorset = 35;

svPath = '/Volumes/Docs/Projects/head_swimming/sim_data';
%svPath = '/Users/mmchenry/Documents/head_swimming/sim_data';
save([svPath filesep 'coord_data'],'d');
save([svPath filesep 'seq_params'],'p');




%% FUNCTIONS

function S = localSystem(P1,P2)
% Defines a rotation matrix for a local coordinate system in an
% inertial frame of reference.  Uses P1 as the origin and P2 to find the
% direction of the x-axis.  Coordinates must be (1x2) vectors.

if size(P1,1)~=1 || size(P1,2)~=2 ||...
   size(P2,1)~=1 || size(P2,2)~=2
    error('Coordinates must be (1x2) vectors');
end

xAxis       = (P2-P1)./norm(P2-P1);
yAxis       = [xAxis(2);-xAxis(1)];
S           = [xAxis' yAxis];


function pts = localToGlobal(pts,origin,S)

if size(pts,2)~=2 || size(origin,2)~=2 
    error('Coordinates must be a (nx2) vector');
end

pts         = [inv(S)'*pts']';
pts(:,1)    = pts(:,1)+origin(1);
pts(:,2)    = pts(:,2)+origin(2);

function pts = globalToLocal(pts,origin,S)

if size(pts,2)~=2 || size(origin,2)~=2 
    error('Coordinates must be a (nx2) vector');
end

pts(:,1)    = pts(:,1)-origin(1);
pts(:,2)    = pts(:,2)-origin(2);
pts         = [S'*pts']';


