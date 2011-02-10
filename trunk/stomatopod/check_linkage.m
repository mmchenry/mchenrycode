function L = check_linkage(p,visData,numVals)
% Evaluates whether the requested geometry is possible and what its limits
% are in range of motion

% Set default values
if nargin < 3
    numVals = 1000;
    if nargin < 2
        visData = 0;
    end
end

% Define full range of possible values for gamma
theta = linspace(0,pi,numVals)';

% Calculate length between B & D for range of theta
h_BD     = sqrt(p.L1^2 + p.L2^2 - 2.*p.L1.*p.L2.*cos(theta));

% Calculate length between B & D for initial theta
h_BD0     = sqrt(p.L1^2 + p.L2^2 - 2.*p.L1.*p.L2.*cos(p.thetaStart));

% Check range of possible gamma
idx  = (h_BD < (p.L3 + p.L4)) & ...
       (h_BD > abs(p.L4-p.L3)) & ...
       (h_BD < (p.L1 + p.L2)) & ...
       (h_BD > abs(p.L1 - p.L2));

% Quit, if no range of motion is possible
if max(idx) == 0  
    warning('No input angles possible with this geometry');
    L = [];
    return
    
elseif ~(  (p.L4 < (p.L3+h_BD0)) || (p.L4 > (h_BD0-p.L3))  )
    warning('Initial linkage geometry not possible');
    L = [];
    return  
end

% Trim impossible range of motion values
theta = theta(idx);
h_BD  = h_BD(idx);

clear idx

% phi - the angle between links 3 and 4
phi  = acos( (p.L3.^2 + p.L4.^2 - h_BD.^2) ./ (2.*p.L3.*p.L4) );

% si - the angle between links 1 and 4
si = acos((h_BD.^2 + p.L1^2 - p.L2^2)./(2*h_BD.*p.L1)) + ...
       acos((h_BD.^2 + p.L4^2 - p.L3^2)./(2*h_BD.*p.L4));

% Positions of points
A = zeros(length(theta),2);
B = [p.L2.*sin(theta) p.L2.*cos(theta)];
C = [p.L4.*sin(si) p.L1-p.L4.*cos(si)];
D = [zeros(length(theta),1) p.L1.*ones(length(theta),1)];

%gamma = pi - atan2((B(:,2)-C(:,2)),-(B(:,1)-C(:,1)));

% Numerical calculation of gamma
gamma = 3*pi/2 + atan2((B(:,1)-C(:,1)),-(B(:,2)-C(:,2)));

% Discrete calculation of KT
cs = fnder(csapi(theta,gamma));
KT = fnval(cs,theta);

% Visualize spline fit used to calculate KT numerically
if 0
    figure
    xVals = theta(2:end)-diff(theta)./2;
    plot(xVals,diff(gamma)./diff(theta),'-k', ...
         xVals,fnval(cs,xVals),'r--')
     clear xVals
end

clear cs

% Find the in vivo values for the gamma range of motion     
gamma_start = gamma(1);   
gamma_end   = gamma(end);   
                         
% Store results
L.theta       = theta;
L.gamma       = gamma;
L.KT          = KT;
L.gamma_start = gamma_start;
L.gamma_end   = gamma_end;


% Visualize curve fit
if visData
    figure;
    plot(L.theta.*(180/pi),L.gamma.*(180/pi),'-',...
         L.linFitIn.*(180/pi),L.linFitOut.*(180/pi),'r-')
    xlabel('theta in');
    ylabel('theta out');
    title(['KT = ' num2str(L.KT_single)]);
    grid on
end




 
     

     
 