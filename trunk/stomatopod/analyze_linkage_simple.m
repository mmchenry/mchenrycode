function analyze_linkage
% Sensitivity analysis that examines KT and Range of Motion (ROM) for 
% different lengths in 4-bar linkage.


% Individual to simulate
indiv = 121;

% Range of length values to consider
dL = 4e-4;

% Load default geometry
p = get_params(indiv);

% Number of length values to consider
numPts = 100;

% Set range of values to examine
L1M = linspace(p.L1-dL,p.L1+dL,numPts);
L2M = linspace(p.L2-dL,p.L2+dL,numPts);
L3M = linspace(p.L3-dL,p.L3+dL,numPts);
L4M = linspace(p.L4-dL,p.L4+dL,numPts);


if 0
l_txt{1} = 'L1';
l_txt{2} = 'L2';
l_txt{3} = 'L3';
l_txt{4} = 'L4';

% Vary L1
j = 1;
for i = 1:length(L1M)
    p.L1 = L1M(i);
    
    L = check_linkage(p,0,10^6);
    
    if isempty(L)
        KTvals(i,j) = nan;
        maxOutVals(i,j) = nan;
    else
        iKT_rest        = find(abs(L.thetaIn-p.thetaRest) == ...
                               min(abs(L.thetaIn-p.thetaRest)));
        %KTvals(i,j)     = L.KT_all(iKT_rest);
        KTvals(i,j)     = min(L.KT_all);
        maxOutVals(i,j) = (L.thetaInMax - L.thetaInMin) * (180/pi);
    end 
    clear L
end

% Vary L2
j = 2;
for i = 1:length(L2M)
    p.L2 = L2M(i);
    
    L = check_linkage(p,0,10^6);
    
    if isempty(L)
        KTvals(i,j) = nan;
        maxOutVals(i,j) = nan;
    else
        iKT_rest        = find(abs(L.thetaIn-p.thetaRest) == ...
                               min(abs(L.thetaIn-p.thetaRest)));
        %KTvals(i,j)     = L.KT_all(iKT_rest);
        KTvals(i,j)     = min(L.KT_all);
        maxOutVals(i,j) = (L.thetaInMax - L.thetaInMin) * (180/pi);
    end 
    clear L
end
 
% Vary L3
j = 3;
for i = 1:length(L3M)
    p.L3 = L3M(i);
    
    L = check_linkage(p,0,10^6);
    
    if isempty(L)
        KTvals(i,j) = nan;
        maxOutVals(i,j) = nan;
    else
        iKT_rest        = find(abs(L.thetaIn-p.thetaRest) == ...
                               min(abs(L.thetaIn-p.thetaRest)));
        %KTvals(i,j)     = L.KT_all(iKT_rest);
        KTvals(i,j)     = min(L.KT_all);
        maxOutVals(i,j) = (L.thetaInMax - L.thetaInMin) * (180/pi);
    end 
    clear L
end

% Vary L4
j = 4;
for i = 1:length(L4M)
    p.L4 = L4M(i);
    
    L = check_linkage(p,0,10^6);
    
    if isempty(L)
        KTvals(i,j) = nan;
        maxOutVals(i,j) = nan;
    else
        iKT_rest        = find(abs(L.thetaIn-p.thetaRest) == ...
                               min(abs(L.thetaIn-p.thetaRest)));
        %KTvals(i,j)     = L.KT_all(iKT_rest);
        KTvals(i,j)     = min(L.KT_all);
        maxOutVals(i,j) = (L.thetaInMax - L.thetaInMin) * (180/pi);
    end 
    clear L
end


% Load default geometry
p = get_params(indiv);

%% Plot

figure;
subplot(2,1,1)
plot(1000.*(L1M-p.L1),KTvals(:,1),'r',...
     1000.*(L2M-p.L2),KTvals(:,2),'b',...
     1000.*(L3M-p.L3),KTvals(:,3),'g',...
     1000.*(L4M-p.L4),KTvals(:,4),'m')
xlabel('Length change (mm)')
ylabel('min KT')
legend(l_txt,'Location','NorthEast')
axis square
xlim([-.5 .5])

subplot(2,1,2)
plot(1000.*(L1M-p.L1),maxOutVals(:,1),'r',...
     1000.*(L2M-p.L2),maxOutVals(:,2),'b',...
     1000.*(L3M-p.L3),maxOutVals(:,3),'g',...
     1000.*(L4M-p.L4),maxOutVals(:,4),'m')
xlabel('Length change (mm)')
ylabel('range of motion')
axis square
xlim([-.5 .5])

end


%% Visualize effect of change in L3

figure;

p = get_params(indiv);
%p.L3 = L3M(1);
L = check_linkage(p,0,10^6);

subplot(2,2,1)
visMorph(p,p.thetaStart,'L3','r')

subplot(2,2,2)
visMorph(p,p.thetaRest,'L3','r')

p.L3 = p.L3+0.0005;
L = check_linkage(p,0,10^6);

subplot(2,2,3)
visMorph(p,p.thetaStart,'L3','r')

subplot(2,2,4)
visMorph(p,p.thetaRest,'L3','r')



%% Visualize 

p = get_params(indiv);
%p.L3 = L3M(1);
L = check_linkage(p,0,10^6);

idx = (L.thetaIn > p.thetaStart) & (L.thetaIn < p.thetaRest);

figure

subplot(2,1,1)
h = plot(L.thetaIn.*(180/pi),L.thetaOut.*(180/pi),'k');
hold on
set(h,'Color',.8.*[1 1 1]);
set(h,'LineWidth',1);
h = plot(L.thetaIn(idx).*(180/pi),L.thetaOut(idx).*(180/pi),'r');

subplot(2,1,2)
h = plot(L.thetaIn.*(180/pi),L.KT_all,'k');
hold on
set(h,'Color',.8.*[1 1 1]);
set(h,'LineWidth',1);
h = plot(L.thetaIn(idx).*(180/pi),L.KT_all(idx),'r');

clear L p

p = get_params(indiv);
p.L3 = p.L3+0.0005;
L = check_linkage(p,0,10^6);

idx = (L.thetaIn > p.thetaStart) & (L.thetaIn < p.thetaRest);

subplot(2,1,1)
h = plot(L.thetaIn.*(180/pi),L.thetaOut.*(180/pi),'k');
hold on
set(h,'Color',.8.*[1 1 1]);
set(h,'LineWidth',1);
h = plot(L.thetaIn(idx).*(180/pi),L.thetaOut(idx).*(180/pi),'b');
set(gca,'YTick',[0:45:180])
axis square

subplot(2,1,2)
h = plot(L.thetaIn.*(180/pi),L.KT_all,'k');
hold on
set(h,'Color',.8.*[1 1 1]);
set(h,'LineWidth',1);
h = plot(L.thetaIn(idx).*(180/pi),L.KT_all(idx),'b');
ylim([2 25])
axis square



function visMorph(p,theta,link,clr)

% Distance between B & D
h_BD  = sqrt(p.L1^2 + p.L2^2 - 2*p.L1*p.L2*cos(theta));

% Angl btwn L4 & L1
si = acos((h_BD.^2+p.L1^2-p.L2^2)/(2*h_BD*p.L1)) + ...
    acos((h_BD^2+p.L4^2-p.L3^2)/(2*h_BD*p.L4));

% Distance between points B & F
h_BF    = sqrt(p.h_AF^2 + p.L2^2 + 2*p.h_AF*p.L2*cos(theta));

% Angle btwn 5 and the x-axis
gamma = pi - acot(p.L5/sqrt(h_BF^2-p.L5^2))...
    - atan(cot(theta) + p.h_AF*csc(theta)/p.L2);

% Define points
A = [0 0];
B = [p.L2*sin(theta) p.L2*cos(theta)];
C = [p.L4*sin(si) p.L1-p.L4*cos(si)];
D = [0 p.L1];
F = [0 -p.h_AF];


% Plot
h = plot([A(1) B(1)],[A(2) B(2)],'k',...
         [B(1) C(1)],[B(2) C(2)],'k',...
         [C(1) D(1)],[C(2) D(2)],'k');
     
set(h,'Color',.8.*[1 1 1]);
set(h,'LineWidth',3);
hold on

if strcmp(link,'L3')
    h = plot([B(1) C(1)],[B(2) C(2)],clr);
else
    h = plot([A(1) B(1)],[A(2) B(2)],clr);
end

set(h,'LineWidth',2);

axis equal
hold off

