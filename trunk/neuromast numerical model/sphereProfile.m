function sphereProfile


%General parameters
mu      = 1.002e-3;
rho     = 998;

f       = 50;
a       = 1.5e-3;


omega   = 2*pi*f;

%Sphere flow parameters
p       = a + 10e-3;
spdSnd  = speedofSound(30);
beta    = (omega .* rho ./(2*mu)).^0.5;
delta   = (2*mu/(rho*omega))^0.5;
h       = (1-i) .* beta;
k       = omega ./ spdSnd;

z       = linspace(p/100000,40e-6,100);

A1      = ((3+3.*i.*h.*a - h.^2.*a.^2) .* k.^3.*a.^3.*exp(i.*k.*a) ./ ...
          (k.^2.*a.^2 .* (1+i.*h.*a) + ...
          (2+2.*i.*k.*a-k.^2.*a.^2) .* h.^2.*a.^2));

Ufull   = lambFlow(a,p-z,h,p,A1,k) - lambFlow(a,p+z,h,p,A1,k);
      
Ulin    = A1 .* exp(-i.*k.*p) .* (-2./(k.*p.^2)+ 6.*i./(k.^2.*p.^3) + ...
          6./(k.^3.*p.^4)) .* z;
      
Upress  = 1 - exp(-z.*(1+i)/delta);

% figure;
% subplot(2,1,1)
% plot(abs(Ufull),z);hold on
% plot(abs(Ulin),z,'r--');
% title([' Lamb (mirror images)  f = ' num2str(f)])
% subplot(2,1,2)
% plot(180*angle(Ufull)/pi,z);hold on
% plot(180*angle(Ulin)/pi,z,'r--')
% 
% disp(' '); disp(' ');
% disp(['tip difference = ' num2str(abs(abs(Ufull(end))-abs(Ulin(end)))/max(z)) ' %'])
% disp(['phase diff = ' num2str(180*(angle(Ufull(end))-angle(Ulin(end)))/pi) ' deg'])
% 
% figure;
% subplot(2,1,1)
% plot(abs(Upress),z,'b--')
% title(['Pressure field  f = ' num2str(f)])
% subplot(2,1,2)
% plot(180*angle(Upress)/pi,z,'b--')



%Normalized versions

Ulin    = Ulin/Ulin(end);
      
Upress  = Upress/Upress(end);


figure;
subplot(2,1,1)
plot(abs(Ulin),z,'r-');hold on
plot(abs(Upress),z,'b--')
subplot(2,1,2)
plot(180*angle(Ulin)/pi,z,'r-'); hold on
plot(180*angle(Upress)/pi,z,'b--')




return
      




figure;
subplot(2,1,1)
plot(abs(U),z);hold on
plot(abs(Upress),z,'b--');
subplot(2,1,2)
plot(180*angle(U)/pi,z);hold on
plot(180*angle(Upress)/pi,z,'b--');


function c = speedofSound(T)
%Gives speed of sound in m/s for temperature in degrees C
% From Bilaniuk & Wong (1993 & 1996)
c =   1.40238744 * 10^3 + 5.03836171 * T - 5.81172916 * 10^-2 * T^2 ...
    + 3.34638117 * 10^-4 * T^3 - 1.48259672 * 10^-6 * T^4 ...
    + 3.16585020 * 10^-9 * T^5;




function u = lambFlow(a,r,h,p,A1,k)

u       = A1 .* (1./(k.^2.*r.^2) + 1./(k.^3.*r.^3)) .* exp(-i.*k.*r);















function c = c_default_theo
%Parameters for all anlayses
 c.freqs         = [10.^linspace(-2,3,100)]';
 c.numHeights    = 50;  
 c.bunHeight     = 5.3e-6; %From Dinklo, 2005
 c.dispAmp       = 10 * 10^-6; %m
 c.E_matrix      = 31; %31 Pa
 c.EI_kino       = 2e-21; % 2e-21 N m^2
 c.bundleStiff   = 2.925e-14; %Nm/rad (van Netten & Kroese, 1987) 
 c.linStiff      = 0.13 * 10^-3; %N/m (van Netten & Kroese, 1987) 
 c.rho           = 998; %998 kg m^-3
 c.mu            = 1.002e-3; %1.002e-3 Pa s

%Data from morphometric measurements (based on stiffness paper)
c.baseDiameter 	= 8.88e-6;
c.midDiameter 	= 7.2e-6 ;
%c.kinoHeight 	= 29.7e-6;
c.kinoHeight 	= 16e-6;
c.cupHeight     = 45e-6;
c.numHairs      = 11;