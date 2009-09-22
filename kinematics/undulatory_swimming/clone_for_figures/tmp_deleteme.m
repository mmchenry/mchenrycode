function tmp_deleteme
% Example data from Chp 17 of Zar
%X = [3 4 5 6 8 9 10 11 12 14 15 16 17]';
%Y = [1.4 1.5 2.2 2.4 3.1 3.2 3.2 3.9 4.1 4.7 4.5 5.2 5.0]';
% 
% sig_x2 = sum(X.^2) - sum(X)^2/length(X);
% sig_y2 = sum(Y.^2) - sum(Y)^2/length(Y);
% sig_xy = sum(X.*Y) - sum(X).*sum(Y)/length(X);
% 
% ss_tot = sig_y2;
% ss_regress = sig_xy^2/sig_x2;
% 
% r2 = ss_regress/ss_tot;



x=[0:.01:1]';y = 5.*x + rand(size(x))-.5;plot(x,y,'.',x,16.*x,'r-')

r2 = rSquared(y,5.*x)


yyy=3;

function r2 = rSquared(X,Y)
% rsquared(X,Y).  This finds the r-squared value that describes the goodness of fit
% for a correlation between column vectors X and Y
if size(X,2)>1 | size(Y,2)>1
	error('You need to use two column vectors!');
end

% Example data from Chp 17 of Zar
%X = [3 4 5 6 8 9 10 11 12 14 15 16 17]';
%Y = [1.4 1.5 2.2 2.4 3.1 3.2 3.2 3.9 4.1 4.7 4.5 5.2 5.0]';

sig_x2 = sum(X.^2) - sum(X)^2/length(X);
sig_y2 = sum(Y.^2) - sum(Y)^2/length(Y);
sig_xy = sum(X.*Y) - sum(X).*sum(Y)/length(X);

ss_tot = sig_y2;
ss_regress = sig_xy^2/sig_x2;

r2 = ss_regress/ss_tot;

[corre,p] = corrcoef([X Y]);
r2=corre(1,2).^2;

