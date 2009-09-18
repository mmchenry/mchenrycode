function m = calcMetrics(morph,calData,param)

%CALCULATE MORPHOMETRICS:
[s,h,w,c,dFin,dFinS,vFin,vFinS]  = organizeData(morph,calData,param.numPts_AP,param.tolerance);
s                                = s-s(1);
[segVol,segArea,V,A]             = calcVolumeArea(s,h,w);
[segfinDArea,totFinAreaD]        = calcFinArea(dFinS,dFin);
[segfinVArea,totFinAreaV]        = calcFinArea(vFinS,vFin);

%STORE DATA:
m.body.s        = s; %body position
m.body.h        = h; %height of meat (dist from dorsal to ventral margins)
m.body.w        = w; %width of meat (dist between left and right margins)
m.body.c        = c; %center of meat in verticle dimension
m.body.V        = segVol; %volume at each body segment
m.body.Vtot     = V; %total body volume
m.body.A        = segArea; %surface area at the segment
m.body.Atot     = A; %total surface area for the body
m.Dfin.s        = dFinS; %body position of dorsal fin segments
m.Dfin.h        = dFin; %height of dorsal fin at segment
m.Dfin.A        = segfinDArea; % area of dorsal fin at segment
m.Dfin.Atot     = totFinAreaD; % total area of dorsal fin at segment
m.Vfin.s        = vFinS; %body position of ventral fin segments
m.Vfin.h        = vFin; %height of ventral fin at segment
m.Vfin.A        = segfinVArea; % area of ventral fin at segment
m.Vfin.Atot     = totFinAreaV; % total area of ventral fin at segment



function [segArea,A] = calcFinArea(s,h)
% This is a simplified version of J. Strother's 'bodyVolume' m-file
% Works only given a static desciption of body volume
% M - mass of the body
% x,y,z - position of the center of mass in the 'morpho' system

if size(s,1) > size(s,2) | size(h,1) > size(h,2)
    error('All inputs must be column vectors');
end
% Calculate trunk and tail step sizes (assumes equal steps)
tStep = mean(diff(s));
if (sum(find(abs(diff(s) - tStep) > .001*tStep)))
    error('trunk length is not equally spaced');
end

% Calculation of volume & area
tArea     = DiscreteInt(tStep,h,'simpson');
A         = tArea(end);
segArea   = [0; diff(tArea)'];


function [s,h,w,c,dFin,sAxialD,vFin,sAxialV] = organizeData(morph,calData,numPts,tolerance);
% Use smoothing spline for each periphery
[xL,yL]     = smoothData(morph.periLeft.x,morph.periLeft.y,numPts,tolerance);
[xR,yR]     = smoothData(morph.periRight.x,morph.periRight.y,numPts,tolerance);
w           = abs(yL-yR).* calData.meanCalConst;
[xD,yD]     = smoothData(morph.periDorsal.x,morph.periDorsal.y,numPts,tolerance);
[xV,yV]     = smoothData(morph.periVentral.x,morph.periVentral.y,numPts,tolerance);
c           = mean([yD;yV],1).* calData.meanCalConst;
h           = abs(yD-yV).* calData.meanCalConst;
s           = xD.* calData.meanCalConst;

% Now, find height of dorsal fins:
dFin.x                      = morph.periDorsalFin.x .* calData.meanCalConst;
dFin.y                      = morph.periDorsalFin.y .* calData.meanCalConst;
finOnBodyD                  = find(s > min(dFin.x));
sAxialD                     = s( finOnBodyD );
hD                          = h( finOnBodyD );
while 1==1
    if ( sAxialD(end)+mean(diff(s)) ) > max(dFin.x), break; end
    sAxialD     = [sAxialD sAxialD(end)+mean(diff(s))]; 
    hD           = [hD 0];
end
[fDorsal.x,fDorsal.y]           = smoothFin(dFin.x,dFin.y,tolerance,sAxialD);
dFin                            = fDorsal.y-hD./2;

% Now, find height of ventral fins:
vFin.x                      = morph.periVentralFin.x .* calData.meanCalConst;
vFin.y                      = morph.periVentralFin.y .* calData.meanCalConst;
finOnBodyV                  = find(s > min(vFin.x));
sAxialV                     = s( finOnBodyV );
hV                          = h( finOnBodyV );
while 1==1
    if ( sAxialV(end)+mean(diff(s)) ) > max(vFin.x), break; end
    sAxialV     = [sAxialV sAxialV(end)+mean(diff(s))]; 
    hV          = [hV 0];
end
[fVentral.x,fVentral.y]         = smoothFin(vFin.x,vFin.y,tolerance,sAxialV);
vFin                            = fVentral.y-hV./2;


function [x,y] = smoothData(x,y,numPts,tolerance)
uniques     = find(~diff(x)==0);
x           = x(uniques);
y           = y(uniques);
sp          = spaps(x, y, tolerance);
x           = min(x):(max(x)-min(x))/(numPts-1):max(x);
y           = fnval(sp,x);

function [xN,yN] = smoothFin(x,y,tolerance,sAxial)
uniques     = find(~diff(x)==0);
x           = x(uniques);
y           = y(uniques);
Ds          = mean(diff(sAxial));
sp          = spaps(x, y, tolerance);
yN          = [fnval(sp,sAxial)];
xN          = sAxial;
if 0,figure;plot(x,y,'.',xN,yN,'g');axis equal;end

function [segVol,segArea,V,A] = calcVolumeArea(s,h,w)
% This is a simplified version of J. Strother's 'bodyVolume' m-file
% Works only given a static desciption of body volume
% M - mass of the body
% x,y,z - position of the center of mass in the 'morpho' system

if size(s,1) > size(s,2) | size(h,1) > size(h,2) | size(w,1) > size(w,2)
    error('All inputs must be column vectors');
end
% Calculate trunk and tail step sizes (assumes equal steps)
tStep = mean(diff(s));
if (sum(find(abs(diff(s) - tStep) > .001*tStep)))
    error('trunk length is not equally spaced');
end

% Calculation of volume & area
tVol      = DiscreteInt(tStep,pi.*h.*w,'simpson');
V         = tVol(end);
segVol    = [0; diff(tVol)'];
tArea     = DiscreteInt(tStep,pi.*(2.*(w.^2+h.^2)-((w-h).^2)./2).^.5,'simpson');
A         = tArea(end);
segArea   = [0; diff(tArea)'];
