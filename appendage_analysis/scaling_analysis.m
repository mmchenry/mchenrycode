function scaling_analysis


%% Path definitions

fPath = '/Volumes/Docs/Projects/Patek_project/Thomas_materials/2010 final/females.csv';
mPath = '/Volumes/Docs/Projects/Patek_project/Thomas_materials/2010 final/males.csv';


%% Import raw data


% Read male and female data
f = csvread(fPath);
m = csvread(mPath);

numParams = 15;
fClr = [1 0 0];
mClr = [0 0 1];

figure;

if size(f,2)~=numParams
    error(' ');
end

params = {'tot length','CL','mass','L1','L2','L3','L4','ang start','ang end',...
          'ang_L2_wire','dist_a_wire','dist_b_e','dist_a_f',...
          'k_spring','D'};
iso_slope = [1 1 3 1 1 1 1 0 0 0 1 1 1 0 5];
      
sclVar   = 2;      

for i = 1:numParams
    subplot(5,3,i);
    
    % Female data
    
    indep = log10(f(:,sclVar));
    dep   = log10(f(:,i));
    
    x = [min(indep) max(indep)];
    
    h = plot(indep,dep,'k.');
    set(h,'Color',fClr)
    hold on
    
    % Run regression, tested vs. isometry
    [pF,b,a,a_pred,r2]= regressTest(indep,dep,iso_slope(i));
    
    h = plot(x,b.*x+a,'k-');
    set(h,'Color',fClr)
    
    h = plot(x,iso_slope(i).*x+a_pred,'k-');
    
    clear indep dep a b p r2
    
   % figure
    
    % Male data    
    indep = log10(m(:,sclVar));
    dep   = log10(m(:,i));
    x     = [min(indep) max(indep)];
    
    h = plot(indep,dep,'k.');
    set(h,'Color',mClr)
    hold on
    
    % Run regression, tested vs. isometry
    [pM,b,a,a_pred,r2]= regressTest(indep,dep,iso_slope(i));
    
    h = plot(x,b.*x+a,'k-');
    set(h,'Color',mClr)
    title(['pF = ' num2str(pF) '   pM = ' num2str(pM)])
    ylabel(['log10(' params{i} ')'])
    xlabel(['log10(' params{sclVar} ')'])
    
    %h = plot(x,iso_slope(i).*x+a_pred,'k-');
    
    %[stats,slope,intercept] = reducedMajorAxis(indep,dep,iso_slope(i),0.05,1)
%     makeLogPlot(f(:,2),f(:,i),pred_slopes(i),'CL',params{i},fClr)
%     hold on
%     makeLogPlot(m(:,2),m(:,i),pred_slopes(i),'CL',params{i},mClr)
    
end
      
return


d.f.tot_length  = f(:,1);
d.f.CL          = f(:,2);
d.f.L1          = f(:,3);
d.f.L2          = f(:,4);
d.f.L3          = f(:,5);
d.f.L4          = f(:,6);
d.f.ang_start   = f(:,7).*(pi/180);
d.f.ang_end     = f(:,8).*(pi/180);
d.f.ang_L2_wire = f(:,9).*(pi/180);
d.f.dist_L2_wire= f(:,10);
d.f.dist_a_wire = f(:,11);
d.f.dist_b_e    = f(:,12);
d.f.dist_a_f    = f(:,13);
d.f.k_spring    = f(:,14);
d.f.D           = f(:,15);


% Male data


d.m.tot_length  = m(:,1);
d.m.CL          = m(:,2);
d.m.L1          = m(:,3);
d.m.L2          = m(:,4);
d.m.L3          = m(:,5);
d.m.L4          = m(:,6);
d.m.ang_start   = m(:,7).*(pi/180);
d.m.ang_end     = m(:,8).*(pi/180);
d.m.ang_L2_wire = m(:,9).*(pi/180);
d.m.dist_L2_wire= m(:,10);
d.m.dist_a_wire = m(:,11);
d.m.dist_b_e    = m(:,12);
d.m.dist_a_f    = m(:,13);
d.m.k_spring    = m(:,14);
d.m.D           = m(:,15);





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

   
