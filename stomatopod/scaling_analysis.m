function scaling_analysis


%% Define data for males and females

% Indvidual numbers
indiv = [3 12 120 121 122 123 125 129 132 133 137];

% Body mass
b_mass = [3.577 4.18 5.473 6.58 5.363 6.497 5.42 5.817 6.147 4.78 3.89];

% Sex (1 for female)
sex = logical([1 1 1 0 1 1 0 0 0 0 1]);

% Carapace length (mm)
CL = [13.9967 14.52 15.457 16.877 15.59 16.24 15.89 16.427 17.3 15.31 13.927];

% Total length (mm)
TL = [53.88 58.6 63.4 66.64 65.44 67.74 63.7 63.82 66 59.0 56.07];

% Starting input angle
theta0 = [78.92 85.21 80 78.76 80 80.35 77.58 78.4 77.42 76.63 78.39];

% Resting input angle
thetaR = [92.18 99.44 92.31 91.76 92.49 92.92 93.41 90.63 ...
                       91.53 89.49 91.57];
                   
% Drag torsion index
D = [0.0769 0.074 0.0756 0.0659 0.073 0.0621 0.0682 0.0741 0.0753 ...
          0.0743 0.0746];
      
% Linear spring stiffness (N/m)               
k_lin = 10.^3 .*[43.777 41.376 27.322 33.649 23.271 72.302 38.244 23.122 ...
                 37.159 35.559 49.321];
             
% Angle between link 2 and wire for applying load during materials testing
lambda = (pi/180).*[26.868 27.258 25.787 28.78 28.125 25.847 29.678 ...
                    28.091 27.934 28.617 25.897];  

% Distance between load and point A:
dist_load = 10.^-3 .*[4.218 4.388 4.821 5.476 4.76 5.152 5.011 5.3895 ...
                      5.8346 5.0327 4.153];

% Torsion spring at mV joint (Nm/rad) 
k_spring = k_lin .* dist_load.^2 .* cos(lambda);

% Linkage coordinates (mm)
Link1 = 10^-3 .* [5.95 6.07 6.536 7.598 6.533 6.885 7.104 7.276 7.798 7.023 5.75];
L2 = 10^-3 .* [3.755 3.903 4.227 4.871 4.1097 4.4355 4.324 4.812 4.9987 4.4675 3.627];
L3 = 10^-3 .* [0.425 0.4456 0.5676 0.595 0.6478 0.614 0.651 0.5832 0.7603 0.5871 0.4925];
L4 = 10^-3 .* [6.676 7.335 7.5588 8.596 7.511 7.988 7.91 8.293 8.718 7.839 6.467];


%% Scaling of mechanical parameters

disp(' ')
disp(' Female CL --------------------------------')
b_iso = 1/3;
%[p,b,a,a_pred,r2]= regressTest(log10(b_mass(sex)),log10(CL(sex)),b_iso);
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(sex)'),log10(CL(sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])

disp(' ')
disp(' Male CL --------------------------------')
b_iso = 1/3;
[p,b,a,a_pred,r2]= regressTest(log10(b_mass(~sex)),log10(CL(~sex)),b_iso);
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(~sex)'),log10(CL(~sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])

disp(' ')
disp(' Female L1 --------------------------------')
b_iso = 1/3;
aa = log10(Link1(sex)');
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(sex)'),aa,b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])

disp(' ')
disp(' Male L1 --------------------------------')
aa = log10(Link1(~sex)');
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(~sex)'),aa,b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])
 
  
disp(' ')
disp(' Female L2 --------------------------------')
b_iso = 1/3;
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(sex)'),log10(L2(sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])

disp(' ')
disp(' Male L2 --------------------------------')
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(~sex)'),log10(L2(~sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])
  

disp(' ')
disp(' Female L3 --------------------------------')
b_iso = 1/3;
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(sex)'),log10(L3(sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])

disp(' ')
disp(' Male L3 --------------------------------')
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(~sex)'),log10(L3(~sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])   
  

disp(' ')
disp(' Female L4 --------------------------------')
b_iso = 1/3;
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(sex)'),log10(L4(sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])

disp(' ')
disp(' Male L4 --------------------------------')
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(~sex)'),log10(L4(~sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])
  
disp(' ')
disp(' Female TL --------------------------------')
b_iso = 1/3;
[p,b,a,a_pred,r2]= regressTest(log10(b_mass(sex)),log10(TL(sex)),b_iso);
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(sex)'),log10(TL(sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])

disp(' ')
disp(' Male TL --------------------------------')
b_iso = 1/3;
[p,b,a,a_pred,r2]= regressTest(log10(b_mass(~sex)),log10(TL(~sex)),b_iso);
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(~sex)'),log10(TL(~sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])
  
disp(' ')
disp(' Female theta0 --------------------------------')
b_iso = 0;
[p,b,a,a_pred,r2]= regressTest(log10(b_mass(sex)),log10(theta0(sex)),b_iso);
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(sex)'),log10(theta0(sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])

disp(' ')
disp(' Male theta0 --------------------------------')
[p,b,a,a_pred,r2]= regressTest(log10(b_mass(~sex)),log10(theta0(~sex)),b_iso);
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(~sex)'),log10(theta0(~sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])
  
disp(' ')
disp(' Female thetaR --------------------------------')
b_iso = 0;
[p,b,a,a_pred,r2]= regressTest(log10(b_mass(sex)),log10(thetaR(sex)),b_iso);
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(sex)'),log10(thetaR(sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])

disp(' ')
disp(' Male thetaR --------------------------------')
[p,b,a,a_pred,r2]= regressTest(log10(b_mass(~sex)),log10(thetaR(~sex)),b_iso);
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(~sex)'),log10(thetaR(~sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])
  
disp(' ')
disp(' Female k_spring --------------------------------')
b_iso = 0;
[p,b,a,a_pred,r2]= regressTest(log10(b_mass(sex)),log10(k_spring(sex)),b_iso);
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(sex)'),log10(k_spring(sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])

disp(' ')
disp(' Male k_spring --------------------------------')
[p,b,a,a_pred,r2]= regressTest(log10(b_mass(~sex)),log10(k_spring(~sex)),b_iso);
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(~sex)'),log10(k_spring(~sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])
  
  
disp(' ')
disp(' Female D --------------------------------')
b_iso = 0;
[p,b,a,a_pred,r2]= regressTest(log10(b_mass(sex)),log10(D(sex)),b_iso);
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(sex)'),log10(D(sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])

disp(' ')
disp(' Male D --------------------------------')
[p,b,a,a_pred,r2]= regressTest(log10(b_mass(~sex)),log10(D(~sex)),b_iso);
[stats,slope,intercept] = ...
       reducedMajorAxis(log10(b_mass(~sex)'),log10(D(~sex)'),b_iso,.05,0);
disp(['a: mean= ' num2str(mean([stats.lowerLimit_a stats.upperLimit_a])) ' L1= ' num2str(stats.lowerLimit_a) ...
      ' L2= ' num2str(stats.upperLimit_a)])
disp(['b: mean= ' num2str(mean([stats.lowerLimit_b stats.upperLimit_b])) ' L1= ' num2str(stats.lowerLimit_b) ...
      ' L2= ' num2str(stats.upperLimit_b)])
  


figure

%% Linear link scaling for females

b_mass = [3.576666667 4.18 5.473333 5.3633 6.49667 3.89];
L1 = [5.86239 6.04197 6.53586 6.53275 6.8850 5.7502];
L2 = [3.755 3.9027 4.2265 4.1097 4.435 3.62734];
L3 = [0.42509 0.44565 0.56764 0.64778 0.61397 0.492494];
L4 = [6.5306 7.29727 7.55877 7.5109 7.988 6.4670];
    
% Run regression, tested vs. isometry
[p,b,a,a_pred,r2]= regressTest(log10(b_mass),log10(L1),1/3);

[stats,slope,intercept] = reducedMajorAxis(log10(b_mass'),log10(L1'),1/3,.05,0);


% Plot L1
%subplot(2,2,1)
h = plot(log10(b_mass),log10(L1),'sg',...
    log10([min(b_mass) max(b_mass)]),...
    log10([min(b_mass) max(b_mass)]).*slope+intercept,'g--');
set(h(1),'MarkerFaceColor','r')
set(h(1),'MarkerSize',6)
axis square
xlabel('Body mass (g)')
ylabel('L1 (mm)')
title(['p = ' num2str(p)])
hold on
xlim([0.5 0.9])
set(gca,'XTick',[.5:.1:.9])

% L2: Run regression, tested vs. isometry
[p,b,a,a_pred,r2]= regressTest(log10(b_mass),log10(L2),1/3);

[stats,slope,intercept] = reducedMajorAxis(log10(b_mass'),log10(L2'),1/3,.05,0);


% Plot L2
%subplot(2,2,2)
h = plot(log10(b_mass),log10(L2),'sr',...
    log10([min(b_mass) max(b_mass)]),...
    log10([min(b_mass) max(b_mass)]).*slope+intercept,'r--');
set(h(1),'MarkerFaceColor','r')
set(h(1),'MarkerSize',6)
axis square
xlabel('Body mass (g)')
ylabel('L2 (mm)')
title(['p = ' num2str(p)])
hold on
xlim([0.5 0.9])
set(gca,'XTick',[.5:.1:.9])


% L3: Run regression, tested vs. isometry
[p,b,a,a_pred,r2]= regressTest(log10(b_mass),log10(L3),1/3);

[stats,slope,intercept] = reducedMajorAxis(log10(b_mass'),log10(L3'),1/3,.05,0);

% Plot L3
%subplot(2,2,3)
h = plot(log10(b_mass),log10(L3),'sk',...
    log10([min(b_mass) max(b_mass)]),...
    log10([min(b_mass) max(b_mass)]).*slope+intercept,'k--');
set(h(1),'MarkerFaceColor','r')
set(h(1),'MarkerSize',6)
axis square
xlabel('Body mass (g)')
ylabel('L3 (mm)')
title(['p = ' num2str(p)])
hold on
xlim([0.5 0.9])
set(gca,'XTick',[.5:.1:.9])


% L4: Run regression, tested vs. isometry
[p,b,a,a_pred,r2]= regressTest(log10(b_mass),log10(L4),1/3);

[stats,slope,intercept] = reducedMajorAxis(log10(b_mass'),log10(L4'),1/3,.05,0);

% Plot L4
%subplot(2,2,4)
h = plot(log10(b_mass),log10(L4),'sb',...
    log10([min(b_mass) max(b_mass)]),...
    log10([min(b_mass) max(b_mass)]).*slope+intercept,'b--');
set(h(1),'MarkerFaceColor','r')
set(h(1),'MarkerSize',6)
axis square
xlabel('Body mass (g)')
ylabel('L4 (mm)')
title(['p = ' num2str(p)])
hold on
xlim([0.5 0.9])
set(gca,'XTick',[.5:.1:.9])


%% Linear link scaling for males

b_mass = [6.58 5.42 5.81667 6.14667 4.78];
L1 = [7.5983 7.104 7.2759 7.798 7.023];
L2 = [4.871 4.3242 4.811 4.99869 4.4675];
L3 = [0.59496 0.650867 0.5832 0.7603 0.587];
L4 = [8.596 7.911 8.293 8.7183 7.8393];
    
% Run regression, tested vs. isometry
[p,b,a,a_pred,r2]= regressTest(log10(b_mass),log10(L1),1/3);

[stats,slope,intercept] = reducedMajorAxis(log10(b_mass'),log10(L1'),1/3,.05,0);

% Plot L1
%subplot(2,2,1)
h = plot(log10(b_mass),log10(L1),'og',...
    log10([min(b_mass) max(b_mass)]),...
    log10([min(b_mass) max(b_mass)]).*slope+intercept,'g-');
set(h(1),'MarkerFaceColor','k')
set(h(1),'MarkerSize',6)
axis square
xlabel('Body mass (g)')
ylabel('L1 (mm)')
title(['p = ' num2str(p)])

% L2: Run regression, tested vs. isometry
[p,b,a,a_pred,r2]= regressTest(log10(b_mass),log10(L2),1/3);

[stats,slope,intercept] = reducedMajorAxis(log10(b_mass'),log10(L2'),1/3,.05,0);

% Plot L2
%subplot(2,2,2)
h = plot(log10(b_mass),log10(L2),'or',...
    log10([min(b_mass) max(b_mass)]),...
    log10([min(b_mass) max(b_mass)]).*slope+intercept,'r-');
set(h(1),'MarkerFaceColor','k')
set(h(1),'MarkerSize',6)
axis square
xlabel('Body mass (g)')
ylabel('L2 (mm)')
title(['p = ' num2str(p)])


% L3: Run regression, tested vs. isometry
[p,b,a,a_pred,r2]= regressTest(log10(b_mass),log10(L3),1/3);

[stats,slope,intercept] = reducedMajorAxis(log10(b_mass'),log10(L3'),1/3,.05,0);

% Plot L3
%subplot(2,2,3)
h = plot(log10(b_mass),log10(L3),'ok',...
    log10([min(b_mass) max(b_mass)]),...
    log10([min(b_mass) max(b_mass)]).*slope+intercept,'k-');
set(h(1),'MarkerFaceColor','k')
set(h(1),'MarkerSize',6)
axis square
xlabel('Body mass (g)')
ylabel('L3 (mm)')
title(['p = ' num2str(p)])


% L4: Run regression, tested vs. isometry
[p,b,a,a_pred,r2]= regressTest(log10(b_mass),log10(L4),1/3);

[stats,slope,intercept] = reducedMajorAxis(log10(b_mass'),log10(L4'),1/3,.05,0);

% Plot L4
%subplot(2,2,4)
h = plot(log10(b_mass),log10(L4),'ob',...
    log10([min(b_mass) max(b_mass)]),...
    log10([min(b_mass) max(b_mass)]).*slope+intercept,'b-');
set(h(1),'MarkerFaceColor','k')
set(h(1),'MarkerSize',6)
axis square
xlabel('Body mass (g)')
ylabel('L4 (mm)')
title(['p = ' num2str(p)])




function makeLogPlot(x,y,prediction,xtext,ytext,clr)

% Scatter plot
h = loglog(x,y,'k.');
hold on;

% Run regression stats
[p,b,a,a_predicted,r2]= regressTest(log10(x),log10(y),prediction);
%[p,b,a,a_predicted]= regressTest2(log10(x),log10(y),prediction);

loglog([min(x) max(x)],[(10.^a).*min(x).^b (10.^a).*max(x).^b],'k');
loglog([min(x) max(x)],[(10.^a_predicted).*min(x).^prediction (10.^a_predicted).*max(x).^prediction],'g');
%loglog([min(x) max(x)],[(10.^a).*min(x).^prediction (10.^a).*max(x).^prediction],'k--');
xlabel(xtext); ylabel(ytext);hold off;
if p>.05 
   title(['a= ' num2str(10.^a) ' b= ' num2str(b) ' p= ' num2str(p) ' r^2= ' num2str(r2) ' , NOT sig. diff. from b=' num2str(prediction)]);
else
   title(['a= ' num2str(10.^a) ' b= ' num2str(b) ' p= ' num2str(p) ' r^2= ' num2str(r2) ' , SIG DIFF from b=' num2str(prediction)]);
end



function [p,slope_fit,intrcpt_fit,intrcpt_pred,r2]= regressTest(Xvalues,Yvalues,slope_pred,intrcpt_pred)
%regressTest addresses the question of whether a regression line found by least-squares
%is significantly different from a known slope.  It runs the least-squares fit on the X and Y values given.
%An entry for kIn constrains the intercept to that value.
%
%See Sokal and Rohlf: pp.465, 471, and the modifications mentioned in pp. 474-475. 
%Generally, a regression is significantly different from the given slope if p < 0.05.

% Find least-squares fit for the slope
n       = length(Xvalues);
DOF     = n-1;
sumX    = sum(Xvalues);
sumY    = sum(Yvalues);
meanX   = mean(Xvalues);
meanY   = mean(Yvalues);
sumX2   = sum([Xvalues-meanX].^2);
sumY2   = sum([Yvalues-meanY].^2);
sumXY   = sum([Yvalues-meanY].*[Xvalues-meanX]);
slope_fit= sumXY/sumX2; 

clear sumX sumY meanX meanY sumY2 

% Define the intercept

if nargin > 3
    % Use intercept provided, adjust DOF
    intrcpt_fit = intrcpt_pred;
    DOF = n;
else
    % Intercepts found by least squares using predicted and fit slopes
    intrcpt_pred = (mean(Yvalues)-slope_pred.*mean(Xvalues));
    
    intrcpt_fit = (mean(Yvalues)-slope_fit.*mean(Xvalues));
end

% Explained variation
exSqr= (sumXY.^2)./sumX2;

% Unexplained variation
unexSqr= sum(([[Yvalues] - ([Xvalues].*slope_fit + intrcpt_fit)]).^2);

if unexSqr~=0,
    
    newYs   = ([Xvalues].*slope_fit + intrcpt_fit);
    s2      = unexSqr./(DOF-1);
    sb      = ( s2./ sumX2 ).^.5;
    tval    = (slope_fit-(slope_pred))./sb;
    p       = tcdf(tval,DOF);
    p       = 2 * min(p,1 - p);
    r2= exSqr./(exSqr+unexSqr);
    
else
    
    p  = 0;
    r2 = 1;
    warning('The trend-line perfectly fits the data and this causes problems in the calculations');
    
end

function [p,b,k,r2]= regressTest_force_int(Xvalues,Yvalues,slope,kIn)
% This version of regressTest addresses the question of whether a regression line found by least-squares
% is significantly different from a known slope when the intercept is provided.  
% It runs the least-squares fit on the X and Y values given.
% An entry for k constrains the intercept to that value.
% See Sokal and Rohlf: pp.465, 471, and the modifications mentioned in pp. 474-475. 
% Generally, a regression is significantly different from the given slope if p < 0.05.

n     = length(Xvalues);
DOF   = n-1;
sumX  = sum(Xvalues);
sumY  = sum(Yvalues);
meanX = mean(Xvalues);
meanY = mean(Yvalues);
sumX2 = sum([Xvalues-meanX].^2);
sumY2 = sum([Yvalues-meanY].^2);
sumXY = sum([Yvalues-meanY].*[Xvalues-meanX]);
b     = sumXY/sumX2; %least-squares solution

if nargin>3
    yIntercept= kIn;
    DOF= n;
else
    %Intecept found by least squares
    yIntercept= (mean(Yvalues)-b.*mean(Xvalues)); 
end

exSqr= (sumXY.^2)./sumX2;
unexSqr= sum(([[Yvalues] - ([Xvalues].*b + yIntercept)]).^2);
if unexSqr==0,
    p=0;
    r2=1;
    k=yIntercept;
else
    newYs= ([Xvalues].*b + yIntercept);
    s2= unexSqr./(DOF-1);
    sb= ( s2./ sumX2 ).^.5;
    tval= (b-(slope))./sb;
    p= tcdf(tval,DOF);
    p= 2 * min(p,1 - p);
    k= yIntercept;
    r2= exSqr./(exSqr+unexSqr);
end
