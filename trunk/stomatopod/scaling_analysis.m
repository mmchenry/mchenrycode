function scaling_analysis


figure

%% scaling for females

b_mass = [3.576666667 4.18 5.473333 5.3633 6.49667 3.89];
L1 = [5.86239 6.04197 6.53586 6.53275 6.8850 5.7502];
L2 = [3.755 3.9027 4.2265 4.1097 4.435 3.62734];
L3 = [0.42509 0.44565 0.56764 0.64778 0.61397 0.492494];
L4 = [6.5306 7.29727 7.55877 7.5109 7.988 6.4670];
    
% Run regression, tested vs. isometry
[p,b,a,a_pred,r2]= regressTest(log10(b_mass),log10(L1),1);

% Plot L1
subplot(2,2,1)
h = plot(log10(b_mass),log10(L1),'or',...
    log10([min(b_mass) max(b_mass)]),...
    log10([min(b_mass) max(b_mass)]).*b+a,'r-');
set(h(1),'MarkerFaceColor','r')
set(h(1),'MarkerSize',4)
axis square
xlabel('Body mass (g)')
ylabel('L1 (mm)')
hold on

% L2: Run regression, tested vs. isometry
[p,b,a,a_pred,r2]= regressTest(log10(b_mass),log10(L2),1);

% Plot L2
subplot(2,2,2)
h = plot(log10(b_mass),log10(L2),'or',...
    log10([min(b_mass) max(b_mass)]),...
    log10([min(b_mass) max(b_mass)]).*b+a,'r-');
set(h(1),'MarkerFaceColor','r')
set(h(1),'MarkerSize',4)
axis square
xlabel('Body mass (g)')
ylabel('L2 (mm)')
hold on


% L3: Run regression, tested vs. isometry
[p,b,a,a_pred,r2]= regressTest(log10(b_mass),log10(L3),1);

% Plot L3
subplot(2,2,3)
h = plot(log10(b_mass),log10(L3),'or',...
    log10([min(b_mass) max(b_mass)]),...
    log10([min(b_mass) max(b_mass)]).*b+a,'r-');
set(h(1),'MarkerFaceColor','r')
set(h(1),'MarkerSize',4)
axis square
xlabel('Body mass (g)')
ylabel('L3 (mm)')
hold on


% L4: Run regression, tested vs. isometry
[p,b,a,a_pred,r2]= regressTest(log10(b_mass),log10(L4),1);

% Plot L4
subplot(2,2,4)
h = plot(log10(b_mass),log10(L4),'or',...
    log10([min(b_mass) max(b_mass)]),...
    log10([min(b_mass) max(b_mass)]).*b+a,'r-');
set(h(1),'MarkerFaceColor','r')
set(h(1),'MarkerSize',4)
axis square
xlabel('Body mass (g)')
ylabel('L4 (mm)')
hold on


%% scaling for males

b_mass = [6.58 5.42 5.81667 6.14667 4.78];
L1 = [7.5983 7.104 7.2759 7.798 7.023];
L2 = [4.871 4.3242 4.811 4.99869 4.4675];
L3 = [0.59496 0.650867 0.5832 0.7603 0.587];
L4 = [8.596 7.911 8.293 8.7183 7.8393];
    
% Run regression, tested vs. isometry
[p,b,a,a_pred,r2]= regressTest(log10(b_mass),log10(L1),1);

% Plot L1
subplot(2,2,1)
h = plot(log10(b_mass),log10(L1),'ok',...
    log10([min(b_mass) max(b_mass)]),...
    log10([min(b_mass) max(b_mass)]).*b+a,'k-');
set(h(1),'MarkerFaceColor','k')
set(h(1),'MarkerSize',4)
axis square
xlabel('Body mass (g)')
ylabel('L1 (mm)')

% L2: Run regression, tested vs. isometry
[p,b,a,a_pred,r2]= regressTest(log10(b_mass),log10(L2),1);

% Plot L2
subplot(2,2,2)
h = plot(log10(b_mass),log10(L2),'ok',...
    log10([min(b_mass) max(b_mass)]),...
    log10([min(b_mass) max(b_mass)]).*b+a,'k-');
set(h(1),'MarkerFaceColor','k')
set(h(1),'MarkerSize',4)
axis square
xlabel('Body mass (g)')
ylabel('L2 (mm)')


% L3: Run regression, tested vs. isometry
[p,b,a,a_pred,r2]= regressTest(log10(b_mass),log10(L3),1);

% Plot L3
subplot(2,2,3)
h = plot(log10(b_mass),log10(L3),'ok',...
    log10([min(b_mass) max(b_mass)]),...
    log10([min(b_mass) max(b_mass)]).*b+a,'k-');
set(h(1),'MarkerFaceColor','k')
set(h(1),'MarkerSize',4)
axis square
xlabel('Body mass (g)')
ylabel('L3 (mm)')


% L4: Run regression, tested vs. isometry
[p,b,a,a_pred,r2]= regressTest(log10(b_mass),log10(L4),1);

% Plot L4
subplot(2,2,4)
h = plot(log10(b_mass),log10(L4),'ok',...
    log10([min(b_mass) max(b_mass)]),...
    log10([min(b_mass) max(b_mass)]).*b+a,'k-');
set(h(1),'MarkerFaceColor','k')
set(h(1),'MarkerSize',4)
axis square
xlabel('Body mass (g)')
ylabel('L4 (mm)')


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
