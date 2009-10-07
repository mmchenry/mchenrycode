function playWithInertia

% %% Test on a cone
% 
% r = 6.7; % Radius 
% L = 20;  % Length
% 
% 
% rho = 4.56;
% 
% s = linspace(0,L,50)';
% a = (r/L).*s;
% b = a;
% c = zeros(length(s),1);
% 
% % This is the result of integration along x
% tmp = 0.25.*a.*b.*pi.*(a.^2 + 4.*(c.^2+s.^2));
% 
% I = rho.*trapz(s,tmp)
% 
% % Analytical solution of I for a cone
% mass = rho*(1/3)*pi*r^2*L;
% (3/5) * mass * (r^2/4 + L^2)
% 
% disp(' '); disp(' ')
% 
% % Integration for COM
% ds = mean(diff(s));
% dA = pi .* a .* b;
% COM_y = rho .* s .* dA;
% trapz(s,COM_y)./ mass
% 
% % Analytical solution
% 0.75*L


% %% Test on a ellpsoid centered about the origin
% 
% disp(' '); disp(' ')
% 
% % Axes that define the axes of the ellipsoid along the x, y, and z axes
% % repectively:
% A = 3.4;
% B = 5.6;
% C = 10;
% 
% % Density
% rho = 4.56;
% 
% % Position along the y-axis
% s  = linspace(-B,B,500);
% 
% % The axes of elliptical cross-sections 
% a0 = real(sqrt(A^2 .* (1 - s.^2./B^2)));
% c0 = real(sqrt(C^2 .* (1 - s.^2./B^2)));
% 
% % z-position of elliptical cross-sections
% z0 = zeros(size(s));
% 
% % This is the result of integration along x
% tmp = 0.25.*a0.*c0.*pi.*(a0.^2 + 4.*(z0.^2+s.^2));
% 
% I = rho.*trapz(s,tmp)
% 
% % Analytical solution 
% mass = rho*(4/3)*pi*A*B*C;
% mass * (A^2+B^2)/5


%% Ellpsoid NOT centered about the origin

disp(' '); disp(' ')

% Axes that define the axes of the ellipsoid along the x, y, and z axes
% repectively:
A = 8.6733e-05;
B = 1.6506e-04;
C = 8.6733e-05;

% Coordinates of ellipsoid center
xC = 0;
yC = -7.5e-04;
zC = -6.2417e-04;

% Density
rho = 1.0272e+03;

% Position along the y-axis
s  = linspace(-B+yC,B+yC,500);

% The axes of elliptical cross-sections 
a0 = real(sqrt(A^2 .* (1 - (s-yC).^2./B^2)));
c0 = real(sqrt(C^2 .* (1 - (s-yC).^2./B^2)));

% z-position of elliptical cross-sections
z0 = zC.*ones(size(s));

% This is the result of integration along x
tmp = 0.25.*a0.*c0.*pi.*(a0.^2 + 4.*(z0.^2+s.^2));

I = rho.*trapz(s,tmp)

% Analytical for very high yC and/or xC (& small zC)
mass = rho*(4/3)*pi*A*B*C;
r    = sqrt(xC^2 + yC^2 + zC^2);
mass * r^2

return


%% Test for a rod

s = linspace(0,100,50)';
a = .5.*ones(length(s)-1,1);
b = .5.*ones(length(s)-1,1);
c = zeros(length(s)-1,1);

Y1 = s(1:end-1);
Y2 = s(2:end);
ds = mean(diff(s));
s = s+ds/2;
s = s(1:end-1);


I = trapz((-1/12) .* a .* b .* pi  .* (3.*a.^2 + 4 .*(3.* c.^2 + Y1.^2 + Y1 .* Y2 + Y2.^2)) .*rho .* (Y1 - Y2))

(rho.*pi*(.5^2)*100)*100^2/3


%% Functions

% function myfun(y,a,b,c)
% 
% x0 = sqrt(a^2 .*(1-y./c^2));
% z0 = sqrt(c^2.*(1-y.^2./b^2).*(1-
% 
% % This is the result of integration along x
% %F = 0.25.*x0.*z0.*pi.*(x0.^2 + 4.*(c.^2+y.^2));
% 
% tmp = 0.25.*a.*b.*pi.*(a.^2 + 4.*(c.^2+s.^2));
% 


