function plot_sphere_tf()
%This is Jim's initial equation for freq. resp. of freq. resp. for vibrating
% sphere (not reconciled with boundayr layer equation


r = 3e-2;
v = 0.894e-3/0.998e3;
c = 1482;
%a = 1e-3;
d = 50e-6;

w = 2*pi*logspace(-2, 3, 100);
a = linspace(0.5e-3, 5e-3, 10);

figure;
ax1 = axes;

figure;
ax2 = axes;

for n = 1:length(a)
   u = vibrosphere(r - d, v, c, a(n), w) - ...
       vibrosphere(r + d, v, c, a(n), w);

   loglog(ax1, w, abs(u));
   hold(ax1, 'on');
   
   semilogx(ax2, w, (180/pi)*angle(u));
   hold(ax2, 'on');
end

hold(ax1, 'on');
loglog(ax1, w, abs(u), 'r');
xlabel(ax1, 'Frequency (Hz)');
ylabel(ax1, 'u_{cupula}/u_{sphere} Amplitude');

hold(ax2, 'on');
loglog(ax2, w, (180/pi)*angle(u), 'r');
xlabel(ax2, 'Frequency (Hz)');
ylabel(ax2, 'u_{cupula}/u_{sphere} Phase Angle');

end


function u = vibrosphere(r, v, c, a, w)

% Calculate some constants
k = w / c;
beta = sqrt(w / (2*v));
h = (1 - i) * beta;

% Calculate leading coefficient
A1 =  ...
       (3 + 3*i*h*a - h.^2*a^2).*k.^3*a^3.*exp(i*k*a) ./ ...
  (k.^2*a^2.*(1+i*h*a) + (2 + 2*i*k*a - k.^2*a^2).*h.^2*a^2);

% Calculate value for u
u = A1.*(i./(k.^2*r^2) + 1./(k.^3*r^3)).*exp(-i*k*r);


end