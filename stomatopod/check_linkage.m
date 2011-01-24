function L = check_linkage(p,visData,numVals)
% Evaluates whether the requested geometry is possible and what its limits
% are in range of motion

if nargin < 3
    numVals = 1000;
    if nargin < 2
        visData = 0;
    end
end

% Define full range of possible values for thetaOut
thetaIn = linspace(0,pi,numVals);

% Calculate length between B & D for range of theta
h_BD     = sqrt(p.L1^2 + p.L2^2 - 2.*p.L1.*p.L2.*cos(thetaIn));

% Calculate length between B & D for initial theta
h_BD0     = sqrt(p.L1^2 + p.L2^2 - 2.*p.L1.*p.L2.*cos(p.thetaStart));

% Check range of possible thetaOut
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
thetaIn = thetaIn(idx);
h_BD     = h_BD(idx);

clear idx

% Calculate possible thetaOut values
thetaOut  = acos( (p.L3.^2 + p.L4.^2 - h_BD.^2) ./ (2.*p.L3.*p.L4) );

% Output range of motion
outROM = max(thetaOut) - min(thetaOut);

% Calculate KT.  Found by evaluating at the first derivative of 
% dInput / dOutput for output = pi/2.
KT_single = (p.L1*p.L2* ...
      sqrt(1 - (p.L1^2+p.L2^2-p.L3^2-p.L4^2)^2/(4*p.L1^2*p.L2^2))) ...
      / (p.L3 * p.L4);

% KT_all = (p.L1*p.L2)/(p.L3*p.L4) .* csc(thetaOut) .* ...
%          sqrt(1-((p.L1^2+p.L2^2-p.L3^2-p.L4^2+2*p.L3*p.L4.*cos(thetaOut)).^2) ...
%          ./ (4*p.L1^2 *p.L2^2));

KT_all = (p.L1*p.L2)/(p.L3*p.L4) .* sin(thetaIn) ./ ...
         sqrt(1-((-p.L1^2-p.L2^2+p.L3^2+p.L4^2+2*p.L1*p.L2.*cos(thetaIn)).^2) ...
         ./ (4*p.L3^2 *p.L4^2));

% Find the in vivo values for the thetaOut range of motion     
thetaOut_start = thetaOut(abs(p.thetaStart-thetaIn) == ...
                             min(abs(p.thetaStart-thetaIn)));   
thetaOut_end = thetaOut(abs(p.thetaRest-thetaIn) == ...
                             min(abs(p.thetaRest-thetaIn)));   
                         
% Store results
L.thetaIn     = thetaIn;
L.thetaOut    = thetaOut;
L.linFitIn    = [thetaIn(1); thetaIn(end)];
L.linFitOut   = L.linFitIn.*KT_single - mean(L.linFitIn.*KT_single) + pi/2;
L.KT_single   = KT_single;
L.KT_all      = KT_all;
L.thetaOutMin = min(thetaOut);
L.thetaOutMax = max(thetaOut);
L.thetaInMin  = min(thetaIn);
L.thetaInMax  = max(thetaIn);
L.thetaOut_start = thetaOut_start;
L.thetaOut_end = thetaOut_end;

% Visualize curve fit
if visData
    figure;
    plot(L.thetaIn.*(180/pi),L.thetaOut.*(180/pi),'-',...
         L.linFitIn.*(180/pi),L.linFitOut.*(180/pi),'r-')
    xlabel('theta in');
    ylabel('theta out');
    title(['KT = ' num2str(L.KT_single)]);
    grid on
end


