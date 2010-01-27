function [threeD,m] = calc3D(morph,calData,param,action)

calK                = calData.meanCalConst;
%renders the body of the fish, given the sructure 'morph'
[s,h,w,c,fin]       = organizeData(morph,param.numPts_AP,param.tolerance);
%s                   = s-s(1);
[segVol,segArea,V,A]= calcVolumeArea(s.*calK,h.*calK,w.*calK);
[Xf,Yf,Zf]          = drawFin_old(fin);
[sFin,segAreaFin,Afin] = calcFinArea(Xf.*calK,Zf.*calK);

switch action
case '3d data'
    [s,h,w,c]           = addMouthCap(s,h,w,c);
    [X,Y,Z]             = drawBody(s,h,w,c,param.numPts_circ);
    %[Xc,Yc,Zc]          = drawMouthCap(s(1),h(1),w(1),c(1),param.numPts_circ);
    
    %Data for visualization:
    threeD.bod.X        = (X-s(1)) .* calK;
    threeD.bod.Y        = Y .* calK;
    threeD.bod.Z        = Z .* calK;
    threeD.fin.X        = (Xf-s(1)) .* calK;
    threeD.fin.Y        = Yf .* calK;
    threeD.fin.Z        = Zf .* calK;
    m                   = [];
    
case 'metrics'
    [streamlining,sNaca,wNaca,hNaca] = streamLineParameter(w,h,s);
  %Data for analysis:
    m.body.length       = max(sFin);
    m.body.s            = s .* calK; %body position
    m.body.h            = h .* calK; %height of meat (dist from dorsal to ventral margins)
    m.body.w            = w .* calK; %width of meat (dist between left and right margins)
    m.body.c            = c .* calK; %center of meat in verticle dimension
    m.body.V            = segVol; %volume at each body segment
    m.body.Vtot         = V; %total body volume
    m.body.streamLining = streamlining;
    m.body.A            = segArea; %surface area at the segment
    m.body.Atot         = A; %total surface area for the body
    m.fin.fHeight       = (max(Zf(:))-min(Zf(:))).*calK;
%     [r,c]               = find(Zf==max(Zf(:)))
%     if r(1)==1, r2=2; elseif r(1)==2,r2=1; elseif r(1)==3,r2=4;elseif r(1)==4,r2=3;end
%     thick1              = abs( Zf(r(1),c(1)) - Zf(r2,c(1)) ).*calK;
%     [r,c]               = find(Zf==min(Zf(:)))
%     if r(1)==1, r2=2; elseif r(1)==2,r2=1; elseif r(1)==3,r2=4;elseif r(1)==4,r2=3;end
%     thick2              = abs( Zf(r(1),c(1)) - Zf(r2,c(1)) ).*calK;
    sDor                = find( (sFin / m.body.length) <.75 );
    m.fin.fHeightNoBody = 2 .* max(abs(max(Zf(:,sDor)) - min(Zf(:,sDor)))).*calK;
    
    %m.fin.fHeightNoBody = thick1 + thick2;
    d                   = find(Xf>.9*max(Xf(:)));
    %m.fin.cFinHeight    = (prctile(Zf(d),95) - prctile(Zf(d),5)).*calK;
    m.fin.cFinHeight    = ( max(Zf(d)) - min(Zf(d)) ).*calK;
    m.fin.s             = sFin; %body position of fin segments
    m.fin.segArea       = segAreaFin; 
    m.Dfin.Atot         = Afin; % total area of dorsal fin at segment
    
  %draw the foil body:
    [s,h,w,c]           = addMouthCap(sNaca,hNaca,wNaca,zeros(size(sNaca)));
    [X,Y,Z]             = drawBody(s,h,w,c,param.numPts_circ);
    %[Xc,Yc,Zc]          = drawMouthCap(s(1),h(1),w(1),c(1),param.numPts_circ);
    
    %Data for visualization:
    threeD.bod.X        = (X-s(1)) .* calK;
    threeD.bod.Y        = Y .* calK;
    threeD.bod.Z        = Z .* calK;
    threeD.fin.X        = [];
    threeD.fin.Y        = [];
    threeD.fin.Z        = [];
end
%disp('Data (Body volume,length):'); 
%m.body.Vtot
%m.body.length 


function [s,segArea,A] = calcFinArea(x,y)
% This is a simplified version of J. Strother's 'bodyVolume' m-file
% Works only given a static desciption of body volume
% M - mass of the body
% x,y,z - position of the center of mass in the 'morpho' system
f           = figure;
h           = fill(x,y,'w');axis equal
yLims       = get(gca,'YLim');
xLims       = get(gca,'XLim');
edges       = findobj(h,'EdgeColor',[0 0 0]);
                  set(edges,'EdgeColor','w')
                  set(gca,'Color','k');
img             = getframe;
img             = img.cdata;
imBW            = im2bw(img,.5);
calCon          = ( max(yLims)-min(yLims) )./size(img,1);
finHeights      = sum(imBW,1).*calCon;
tStep           = calCon;
s               = (1:size(img,2)).*calCon + min(xLims);

% Calculation of area
tArea     = DiscreteInt(tStep,finHeights,'simpson');
A         = tArea(end);
segArea   = [0; diff(tArea)'];

%Account for 2 sides of fin:
segArea   = 2 .* segArea;
if 0,hold on,plot(s,finHeights,'m');end
close(f);

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


function [s,h,w,c] = addMouthCap(s,h,w,c)
Ds      = (s(2)-s(1))./10;
s       = [s(1)-2*Ds s(1)-Ds  s];
h       = [h(1)./100 h(1)./2 h];
w       = [w(1)./100 w(1)/2 w];
c       = [c(1) c(1) c];

function [s,h,w,c,fin] = organizeData(morph,numPts,tolerance);
% Use smoothing spline for each periphery
    [xL,yL]     = smoothData(morph.periLeft.x,morph.periLeft.y,numPts,tolerance);
    [xR,yR]     = smoothData(morph.periRight.x,morph.periRight.y,numPts,tolerance);
    w           = abs(yL-yR);
    [xD,yD]     = smoothData(morph.periDorsal.x,morph.periDorsal.y,numPts,tolerance);
    [xV,yV]     = smoothData(morph.periVentral.x,morph.periVentral.y,numPts,tolerance);
    c           = mean([yD;yV],1);
    h           = abs(yD-yV);
    s           = xD;
% Now, deal with the fins
    finTemp.x               = [morph.periDorsalFin.x; morph.periVentralFin.x(length(morph.periVentralFin.x):-1:1)];
    finTemp.y               = [morph.periDorsalFin.y; morph.periVentralFin.y(length(morph.periVentralFin.y):-1:1)];
    [fin.tip.x,fin.tip.y]   = smoothFin(finTemp.x,finTemp.y,2*numPts,tolerance);
    fin.tip.x               = fin.tip.x;
    fin.tip.y               = fin.tip.y;
    bodTemp.x               = [morph.periDorsal.x; morph.periVentral.x(length(morph.periVentral.x):-1:1)];
    bodTemp.y               = [morph.periDorsal.y; morph.periVentral.y(length(morph.periVentral.y):-1:1)];
    finOnBody               = find(bodTemp.x>min(fin.tip.x));
    bodTemp.x               = bodTemp.x(finOnBody);
    bodTemp.y               = bodTemp.y(finOnBody);
    [fin.base.x,fin.base.y] = smoothFin(bodTemp.x,bodTemp.y,2*numPts,tolerance);
    %fin.base.x              = fin.base.x - s(1);
    %fin.tip.x               = fin.base.x - s(1);
    %s                       = s - s(1);
    
    %figure;fill(finData.x,finData.y,[1 0 0]);axis equal

function [x,y] = smoothData(x,y,numPts,tolerance)
uniques     = find(~diff(x)==0);
x           = x(uniques);
y           = y(uniques);
sp          = spaps(x, y, tolerance);
x           = min(x):(max(x)-min(x))/(numPts-1):max(x);
y           = fnval(sp,x);

function [x,y] = smoothFin(x,y,numPts,tolerance)
uniques     = find(~diff(x)==0);
%x           = x(uniques);
%y           = y(uniques);
s           = [0; cumsum( (diff(x).^2 + diff(y).^2).^.5 )];
sNew        = 0:max(s)./(numPts-1):max(s);
m           = [x y];
sp          = spaps(s, m', tolerance);
p           = [fnval(sp,sNew)]';
x           = p(:,1);
y           = p(:,2);
if 0,figure;plot(m(:,1),m(:,2),'.',p(:,1),p(:,2),'g',m(1,1),m(1,2),'o' );axis equal;end

function [x,y,z] = drawFin(fin)
revv    = length(fin.tip.x):-1:1;
x       = [fin.base.x; fin.tip.x(revv)];
z       = [fin.base.y; fin.tip.y(revv)];
y       = zeros(size(x));

function [x,y,z] = drawFin_old(fin)
all     = 1:length(fin.base.x)-1;
x       = [fin.base.x(all)'; fin.tip.x(all)';fin.tip.x(all+1)';fin.base.x(all+1)'];
z       = [fin.base.y(all)'; fin.tip.y(all)';fin.tip.y(all+1)';fin.base.y(all+1)'];
y       = zeros(size(x));

function [X,Y,Z] = drawMouthCap(s,h,w,c,numPts)
%[Ytemp,Ztemp]   = ellipse1(0,0,[1 axes2ecc(w/2,h/2)],0,[0 360],[],'degrees',numPts);
[Ytemp,Ztemp]   = giveEllipse(0,0,w/2,h/2,numPts);
Y               = w/2 * Ytemp;
Z               = w/2 * Ztemp + c;
X               = s*ones(size(Y));

function [x,y,z]= drawBody(s,h,w,c,numPts)
x=[];y=[];z=[];
for i=1:length(s)-1  
  %DRAW FIRST ELLIPSE:
    height          = h(i);
    width           = w(i);
    ecc             = ( 1 - (height/2)^2/(width/2)^2 )^.5;
    [yTemp1,zTemp1] = ellipse1(0,0,[1 ecc],0,[0 360],[],'degrees',numPts);
  % [yTemp1,zTemp1] = giveEllipse(0,c(i),width/2,height/2,numPts);
   % [yTemp1,zTemp1] = ellipse1(0,0,[1 axes2ecc(width/2,height/2)],0,[0 360],[],'degrees',numPts);
    zTemp1          = (width/2) * zTemp1 + c(i);
    yTemp1          = (width/2) * yTemp1;
    xTemp1          = s(i)*ones(size(yTemp1));
  %DRAW SECOND ELLIPSE:
    height          = h(i+1);
    width           = w(i+1);
    centerr         = c(i+1);
  %  [yTemp2,zTemp2] = giveEllipse(0,c(i+1),width/2,height/2,numPts);
    ecc             = ( 1 - (height/2)^2/(width/2)^2 )^.5;
    [yTemp2,zTemp2] = ellipse1(0,0,[1 ecc],0,[0 360],[],'degrees',numPts);
    %[yTemp2,zTemp2] = ellipse1(0,0,[1 axes2ecc(width/2,height/2)],0,[0 360],[],'degrees',numPts);
    zTemp2          = (width/2) * zTemp2 + c(i+1);
    yTemp2          = (width/2) * yTemp2;
    xTemp2          = s(i+1)*ones(size(yTemp2));   
    if 0,plot3(xTemp1,yTemp1,zTemp1,'g',xTemp2,yTemp2,zTemp2,'m');axis equal;hold on;end
  %COMBINE DATA:
    all             = 1:length(xTemp2)-1;
    x               = [x [xTemp1(all)'; xTemp2(all)'; xTemp2(all+1)'; xTemp1(all+1)']];
    y               = [y [yTemp1(all)'; yTemp2(all)'; yTemp2(all+1)'; yTemp1(all+1)']];
    z               = [z [zTemp1(all)'; zTemp2(all)'; zTemp2(all+1)'; zTemp1(all+1)']];
end  

function [x,y] = giveEllipse(xCenter,yCenter,a,b,numPts)
theta       = linspace(0,2*pi,numPts);
x           = a .* cos(theta) + xCenter;
y           = b .* sin(theta) + yCenter;

function makeImage(X,Y,Z,alphaNum,color)
%makes the image
    if isstruct(X)
        p=patch(X);
    else
        p= patch(real(X),real(Y),real(Z),real(Z)*0);
    end
    set(p,'FaceLighting','gouraud','LineStyle','none','BackFaceLighting','reverselit');
    alpha(p,alphaNum);
    set(p,'FaceColor',color);
    axis equal;
    light;
    view([40 25]);
    set(gcf,'Color',[0 0 0]);
    f2= get(gcf,'Children');
    set(f2,'Color',[0 0 0]);
    set(f2,'ZDir','reverse');
    %set(gcf,'Renderer','OpenGL');
    %set(p,'EraseMode','xor');

function viewSide(action)
switch action
case 'dorsal', view([0 90]); case 'ventral', view([0 270]); case 'right', view([0 360]);
case 'left', view([180 360]); case 'headon', view([90 0]); case 'tailon', view([270 0]);
case '3/4', view(3); case '3/4 front',view([40 25]);
end

function fig(action,f,w,h)
if nargin<1, action= 'light black'; end
switch action
case 'light white'
    set(f,'Color',[1 1 1]);
    set(f,'Renderer','OpenGL');
    light;
case 'light black'
    set(f,'Color',[0 0 0]);
    f2= get(f,'Children');
    set(f2,'Color',[0 0 0]);
    set(f,'Renderer','OpenGL');
    light; %set(f2,'nextplot','replace','Visible','off');
case 'flat white'
    set(f,'Color',[1 1 1]);
    f2= get(f,'Children');
    set(f2,'Color',[1 1 1]);
    set(f,'Renderer','OpenGL');  
end
if w==1000 & h==1000
    set(f,'Position',[1030 -250 1000 1000]);%use with pixel dimensions of 1000
else
    set(f,'Position',[1030 1 w h]);%388 408 %543 571
end
set(f,'DoubleBuffer','on');
axis equal;

 function s= localSystem3(p1,p2,p3)
 %defines a local coordinate system (unit vectors) for two points given,p1 is the xaxis, p2 is the origin, p3 is the zaxis
    p1= p1-p2;
    xaxis= p1./norm([p1]);
    zaxis= p3./norm([p3]);
    yaxis= cross(zaxis,xaxis);   
    s= [xaxis' yaxis' zaxis'];