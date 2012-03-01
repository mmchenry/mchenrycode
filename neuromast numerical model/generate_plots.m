function generate_plots

%This was written by Jim

% Physical constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c = 1450; % speed of sound in water (m/s)
c  = 1497; % speed of sound in water at 25C (m/s)
v = 1e-6; % kinematic viscosity of water (m^2/s)
a = 3e-3; % radius of the oscillating sphere (m)
r = 10e-3; % distance from sphere center (m)

% Construct plot showing relative contribution of each term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = logspace(-2, 0, 100); % frequencies (1/s)
%f       = [10.^linspace(-4,0,100)];
[u, ua, ub] = oscillating_sphere(c, v, 2*pi*f, a, r);

figure;
loglog(f, abs(ua), f, abs(ub));
legend('e^{-i k a}', 'e^{-i h a}', 'Location', 'SouthWest');
xlabel('Frequency (Hz)');
ylabel('Relative Ressponse (unitless)');


% Construct plot showing relative contribution of each term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = logspace(-4, 3, 100); % frequencies (1/s)
[u, ua, ub] = oscillating_sphere(c, v, 2*pi*f, a, r);

figure;
subplot(2,1,1);
loglog(f, abs(u));
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot(2,1,2);
semilogx(f, 180*unwrap(angle(u))/pi);
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');


function [u, ua, ub] = oscillating_sphere(c, v, w, a, r)

   k = sqrt(w.^2 ./ (c.^2 + (4/3)*i*v.*w));
   h = (1-i) * sqrt(w ./ (2*v));

   A1 = - (3 + 3*i*h.*a - h.^2.*a.^2).*k.^3.*a.^3.*exp(i*k.*a) ./...
     (k.^2.*a.^2.*(1+i*h.*a) + (2+2*i*k.*a-k.^2.*a.^2).*h.^2.*a.^2);

   B1 = (1/3) * (3 + 3*i*k.*a - k.^2.*a.^2).*h.^3.*a.^3.*exp(i*h.*a) ./...
     (k.^2.*a.^2.*(1+i*h.*a) + (2+2*i*k.*a-k.^2.*a.^2).*h.^2.*a.^2);

   ua = - A1.*(i./(k.^2.*r.^2) + 1./(k.^3.*r.^3)).*exp(-i*k.*r);
   ub = B1.*(3./(h.*r) - 3*i./(h.^2.*r.^2) - 3./(h.^3.*r.^3)).*exp(-i*h.*r);

   u = ua + ub;
