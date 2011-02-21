function [stats,slope,intercept] = reducedMajorAxis(X,Y,b_predicted,alpha,plotting)
% Runs appropriate type 2 regression stats for scaling data according to Rayner, 1985
%plotting=1 makes plots, 0 doesn't
% Code tested using examples from Sokal and Rolf (1995) p.546, where
% x = [14;17;24;25;27;33;34;37;40;41;42];y=[61;37;65;69;54;93;87;89;100;90;97]; 
if nargin < 4, alpha   = .05;end
if nargin < 5, plotting = 1;end

if (size(X,2)>size(X,1)) | (size(Y,2)>size(Y,1))
    error('Data should be arranged in columns, not rows');
end

%In the case that you do not have information about the error of your samples, it is appropriate
%to use a Reduced Major Axis (i.e. geometric mean) regression.  This is the case if you don't have repeated
%measures (i.e. the x and y data have only one column).
if size(X,2)==1
  %First, find standard linear regression stats:
    x       = X - mean(X);
    y       = Y - mean(Y);
    n       = size(Y,1);
    b       = sum(x.*y)./sum(x.^2); %slope
    a       = mean(Y) - b * mean(X); %intercept
    Ypred   = b.*X + a; %predicted y values
    s_b     = ( (sum((Y-Ypred).^2) ./ (n-2)) ./ ...
              sum(x.^2) ).^.5; %standard error of regression
    s_Y     = ( (sum((Y-Ypred).^2) ./ (n-2)) ./ ...
               n ).^.5; %standard error of regression
   %calculate confidence intervals for linear regression:
    tValue  = tinv(1-alpha/2,n-2);
    t_s     = tValue.*s_b;
    t_sa    = tValue.*s_Y;
    L1      = b - t_s;
    L2      = b + t_s;
  %Now, find RMA stats:    
    v       = ( sum(y.^2) / sum(x.^2) ) .^.5;
  %note: the stat v normally cannot find the sign of the slope, 
  %so we use the least squares regression to find the sign:
    c       = corrcoef(X,Y);
    slope   = (c(2)/abs(c(2))) * v;
    a_v     = mean(Y) - slope * mean(X);
    s_v     = s_b;
  %and find confidence intervals for v:
    L1      = v - t_s;
    L2      = v + t_s;
  %and find confidence intervals for A_v:
    aL1     = a_v - t_sa;
    aL2     = a_v + t_sa;
  %prediction stats:
    a_pred  = mean(Y) - b_predicted * mean(X);
    if (b_predicted<L2) && (b_predicted>L1)
        stats.hypothesis    = 'accept null';
    else 
        stats.hypothesis    = 'reject null';
    end
    intercept           = a_v;
    stats.lowerLimit_b  = L1;
    stats.upperLimit_b  = L2;
    stats.lowerLimit_a  = aL1;
    stats.upperLimit_a  = aL2;
    stats.alpha         = alpha;
    if plotting
        %figure;
        pointColor  = 'k';
        h = plot(X,Y,'o');hold on
            set(h,'MarkerFaceColor',pointColor)
            set(h,'MarkerEdgeColor',pointColor)
            set(h,'MarkerSize',3)
        h = plot([min(X) max(X)],slope*[min(X) max(X)]+a_v,'k');
        h = plot([min(X) max(X)],b_predicted*[min(X) max(X)]+a_pred,'g');
        title(stats.hypothesis)
    end
    
    
end