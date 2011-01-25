function define_morph(p,thetaIn)
% Provides coordinates for all points in the morphology for required input
% parameters


%% Error check geometry
   
% Calculate length between B & D for range of theta
h_BD     = sqrt(p.L1^2 + p.L2^2 - 2.*p.L1.*p.L2.*cos(thetaIn));

% Check that linkage geometry is possible
if h_BD > (p.L3 + p.L4)
    error('h_BD cannot exceed L3 + L4')
elseif h_BD < abs(p.L4-p.L3)
    error('h_BD cannot be less than L4 - L3')       
elseif h_BD > (p.L1 + p.L2)
    error('h_BD cannot exceed L1 + L2')
elseif h_BD < abs(p.L1 - p.L2)
    error('h_BD cannot be less than L1 - L2')
end

% Angle between links 1 and 4   
si = acos((h_BD^2 + p.L1^2 - p.L2^2)/(2*h_BD*p.L1)) + 
     acos((h_BD^2 + p.L4^2 - p.L3^2)/(2*h_BD*p.L4)];


(* Define appendage coordinates ------------------------------------------*)
AX = 0;
AY = 0;
BX = L2 Sin[thetaIn];
BY = L2 Cos[thetaIn];
CX = L4 Sin[si];
CY = L1 - L4 Cos[si];
DX = 0;
DY = L1;
FX = 0;
FY = -hAF;
GX = FX - L5 Sin[gamma];
GY = FY + L5 Cos[gamma];

