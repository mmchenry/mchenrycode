function analyze_linkage



%% Prompt for data file, if not given


[fileName,pathh] = uigetfile('_linkdata.mat','Choose a data file');


[fName,fPath,fExt] = giveFileParts(fileName,pathh);

clear fileName fExt


%% Define points in phots's frame of reference (FoR), g
load(fPath)

Ag = [coords.x(1) coords.y(1)];
Bg = [coords.x(2) coords.y(2)];
Cg = [coords.x(3) coords.y(3)];
Dg = [coords.x(4) coords.y(4)];
Eg = [coords.x(5) coords.y(5)];
Fg = [coords.x(6) coords.y(6)];
Gg = [coords.x(7) coords.y(7)];
Hg = [coords.x(8) coords.y(8)];

clear coords


%% Define new point at intersection of AD and CH
% This will be used as the new L4 b/c it is closer to the line of action
% for the muscle

% Coefficients for AD
cAD(1) = (Dg(2)-Ag(2))/(Dg(1)-Ag(1));
cAD(2) = Dg(2)-cAD(1)*Dg(1);

% Coefficients for CH
cCH(1) = (Hg(2)-Cg(2))/(Hg(1)-Cg(1));
cCH(2) = Hg(2)-cCH(1)*Hg(1);

% Solve for coincident point, I
Ig(1,1) = (cCH(2)-cAD(2)) / (cAD(1)-cCH(1));
Ig(1,2) = Ig(1,1) * cCH(1) + cCH(2);

clear cAD cCH


%% Transform all points to merus FoR, m

S_m = localSystem(Ag,Ig);

Am = globalToLocal(Ag,Ag,S_m);
Bm = globalToLocal(Bg,Ag,S_m);
Cm = globalToLocal(Cg,Ag,S_m);
Dm = globalToLocal(Dg,Ag,S_m);
Em = globalToLocal(Eg,Ag,S_m);
Fm = globalToLocal(Fg,Ag,S_m);
Gm = globalToLocal(Gg,Ag,S_m);
Hm = globalToLocal(Hg,Ag,S_m);
Im = globalToLocal(Ig,Ag,S_m);

clear S_m

%% Visually confirm transformation into m FoR

figure;

plot([Ag(1) Bg(1) Cg(1) Dg(1) Ag(1)],[Ag(2) Bg(2) Cg(2) Dg(2) Ag(2)],'m', ...
     [Eg(1) Fg(1) Gg(1) Eg(1)],[Eg(2) Fg(2) Gg(2) Eg(2)],'m',...
     [Bg(1) Eg(1) Cg(1) Bg(1)],[Bg(2) Eg(2) Cg(2) Bg(2)],'m',...
     [Cg(1) Hg(1)],[Cg(2) Hg(2)],'m--',Ag(1),Ag(2),'mo');  hold on
 
plot([Am(1) Bm(1) Cm(1) Dm(1) Am(1)],[Am(2) Bm(2) Cm(2) Dm(2) Am(2)],'b', ...
     [Em(1) Fm(1) Gm(1) Em(1)],[Em(2) Fm(2) Gm(2) Em(2)],'b',...
     [Bm(1) Em(1) Cm(1) Bm(1)],[Bm(2) Em(2) Cm(2) Bm(2)],'b',...
     [Cm(1) Hm(1)],[Cm(2) Hm(2)],'b--',Am(1),Am(2),'bo'); 
 
axis equal 
grid on

% Clear coordinates in photo's frame of reference
clear Ag Bg Cg Dg Eg Fg Gg Hg Ig


%% Define curpus/dactyl points in carpus FoR

S_c = localSystem(Bm,Cm);

Bc = globalToLocal(Bm,Bm,S_c);
Cc = globalToLocal(Cm,Bm,S_c);
Ec = globalToLocal(Em,Bm,S_c)
Fc = globalToLocal(Fm,Bm,S_c)
Gc = globalToLocal(Gm,Bm,S_c)

% Angular position of the carpus in the merus FoR
angl = pi/2 - atan((Cm(2)-Bm(2))/(Cm(1)-Bm(1)))






%% Functions

function S = localSystem(P1,P2)
% Defines a local coordinate system from points in a global frame of
% reference
% P1 (1 x 2) is the origin
% P2 (1 x 2) is the y-axis

yAxis   = (P2-P1)./norm(P2-P1);
xAxis   = [yAxis(2);-yAxis(1)];
S       = [xAxis yAxis'];


function pts = localToGlobal(pts,origin,S)
% Transforms pts (n x 2) from local to global coodinate system

pts     = [inv(S)'*pts']';
pts(:,1)= pts(:,1)+origin(1);
pts(:,2)= pts(:,2)+origin(2);


function pts = globalToLocal(pts,origin,S)
% Transforms pts (n x 2) from global to local coordinate system

pts(:,1)= pts(:,1)-origin(1);
pts(:,2)= pts(:,2)-origin(2);
pts     = [S'*pts']';


function [fileName,filePath,fileExt] = giveFileParts(fileName,pathh)

tmp = find(fileName=='.');

if length(tmp)>1
    error('You cannot have periods in the filename');
elseif length(tmp)<1
    error('Your filename needs an extension')
end

fileExt  = fileName(tmp:end);
fileName = fileName(1:tmp-1);
filePath = [pathh filesep fileName fileExt];