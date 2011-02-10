function [sp_gamma,sp_KT,sp_theta,min_theta,max_theta] = linkage_functions(L1,L2,L3,L4)
% Returns spline functions that describe how gamma and KT vary with theta

% Number of discrete values to consider
numVals = 1000;

% Define full range of possible values for gamma
theta = linspace(0,pi,numVals)';

% Calculate length between B & D for range of theta
h_BD     = sqrt(L1^2 + L2^2 - 2.*L1.*L2.*cos(theta));

% Check range of possible gamma
idx  = (h_BD < (L3 + L4)) & ...
       (h_BD > abs(L4-L3)) & ...
       (h_BD < (L1 + L2)) & ...
       (h_BD > abs(L1 - L2));

% Quit, if no range of motion is possible
if max(idx) == 0  
    error('No input angles possible with this geometry');
end

% Trim impossible range of motion values
theta = theta(idx);
h_BD  = h_BD(idx);

clear idx

% phi - the angle between links 3 and 4
phi  = acos( (L3.^2 + L4.^2 - h_BD.^2) ./ (2.*L3.*L4) );

% si - the angle between links 1 and 4
si = acos((h_BD.^2 + L1^2 - L2^2)./(2*h_BD.*L1)) + ...
       acos((h_BD.^2 + L4^2 - L3^2)./(2*h_BD.*L4));

% Positions of points
A = zeros(length(theta),2);
B = [L2.*sin(theta) L2.*cos(theta)];
C = [L4.*sin(si) L1-L4.*cos(si)];
D = [zeros(length(theta),1) L1.*ones(length(theta),1)];

%gamma = pi - atan2((B(:,2)-C(:,2)),-(B(:,1)-C(:,1)));

% Numerical calculation of gamma
gamma = 3*pi/2 + atan2((B(:,1)-C(:,1)),-(B(:,2)-C(:,2)));

% Spline fit to gamma values
sp_gamma = csapi(theta,gamma);

% Spline fit to theta values
sp_theta = csapi(gamma,theta);

% Find KT, the first derivative of gamma spline
sp_KT = fnder(sp_gamma);

% Range of possible theta values
min_theta = min(theta);
max_theta = max(theta);