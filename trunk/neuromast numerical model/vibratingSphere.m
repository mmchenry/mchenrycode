function vibratingSphere

% Load model parameters and run simulation
c       = c_default_theo;
c       = numerical_twopart(c,'torsion spring');
S_bun   = [c.M(1,:) ./ (i.*(c.freqs').*2.*pi.*c.dispAmp)];

% Extract parameter values
%f       = c.freqs';
f = [10.^linspace(-4,2,100)];
rho     = c.rho;
mu      = c.mu;
zeta    = c.cupHeight;

clear c

% Boundry layer parameters
omega   = f .* 2*pi;
delta   = ( 2 .* mu ./ rho ./ omega ).^0.5;

% Sphere flow parameters
r       = 10e-3; % distance from center of sphere (m)
spdSnd  = speedofSound(20);
beta    = (omega .* rho ./(2*mu)).^0.5;
h       = (1-i) .* beta;
k       = omega ./ spdSnd;

%Check assumptions
deltaMin    = ( 2 .* mu ./ rho ./ max(omega) ).^0.5;

disp(' ');disp(' ');
disp(['zeta/delta = ' num2str(zeta/deltaMin) ' -- Should be << 1']);disp(' ')
disp(['zeta/r = ' num2str(zeta/r) ' -- Should be << 1']);disp(' ')
disp(['zeta*k = ' num2str(zeta*max(k)) ' -- Should be << 1']);disp(' ')



% Plot for range of sphere radii - FLOW DEFINED at bundle height___________
%a  = [1:2:7].*10^-3;
bh =  5.3e-6;
a  = 3.*10^-3;
f1 = figure;

for j = 1:length(a)    
    A1      = ((3+3.*i.*h.*a(j) - h.^2.*a(j).^2) .* k.^3.*a(j).^3.*exp(i.*k.*a(j)) ./ ...
                (k.^2.*a(j).^2 .* (1+i.*h.*a(j)) + ...
                (2+2.*i.*k.*a(j)-k.^2.*a(j).^2) .* h.^2.*a(j).^2));

    U       = A1 .* exp(-i.*k.*r) .* (-2./(k.*r.^2)+ 6.*i./(k.^2.*r.^3) + ...
              6./(k.^3.*r.^4)) .* bh;


    figure(f1)
    subplot(2,1,1)
    loglog(f,abs(U),'b-');hold on
    subplot(2,1,2)
    semilogx(f,angle(U)*180/pi,'b-'); hold on
    pause(.5)

end




% Plot for range of sphere radii_________________________________________
%a  = [1:2:7].*10^-3;
a  = 3.*10^-3;
f1 = figure;
f2 = figure;

disp(['Distance between edge of sphere and surface = ' num2str(r-a)]);
for j = 1:length(a)    
    A1      = ((3+3.*i.*h.*a(j) - h.^2.*a(j).^2) .* k.^3.*a(j).^3.*exp(i.*k.*a(j)) ./ ...
                (k.^2.*a(j).^2 .* (1+i.*h.*a(j)) + ...
                (2+2.*i.*k.*a(j)-k.^2.*a(j).^2) .* h.^2.*a(j).^2));

    E_inf   = (1-i).*A1.*delta .* exp(-i.*k.*r) .* ...
                (-(1./k.*r.^2) + 3*i./(k.^2.*r.^3) + 3./(k.^3.*r.^4));

    
    S_tot   = S_bun .* E_inf;
   

    figure(f1)
    subplot(2,1,1)
    loglog(f,abs(E_inf)); hold on
    subplot(2,1,2)
    semilogx(f,angle(E_inf)*180/pi); hold on


    figure(f2)
    subplot(2,1,1)
    loglog(f,abs(S_tot),'r-');hold on
    subplot(2,1,2)
    semilogx(f,angle(S_tot)*180/pi,'r-'); hold on
    pause(.5)

end

figure(f1)
subplot(2,1,1)
ylabel('Sensitivity (s)')
xlabel('Freq. (Hz)');
modPlot
subplot(2,1,2)
ylabel('Phase (deg)');
xlabel('Freq. (Hz)');
modPlot
    
figure(f2)
subplot(2,1,1)
ylabel('Sensitivity (s)')
xlabel('Freq. (Hz)');
modPlot
subplot(2,1,2) 
ylabel('Phase (deg)');
xlabel('Freq. (Hz)');
modPlot


return


% %Plot for all layers of filtering_________________________________________
% % Plot for range of sphere radii
% a  = 5.*10^-3;
%    
% A1      = ((3+3.*i.*h.*a - h.^2.*a.^2) .* k.^3.*a.^3.*exp(i.*k.*a) ./ ...
%     (k.^2.*a.^2 .* (1+i.*h.*a) + ...
%     (2+2.*i.*k.*a-k.^2.*a.^2) .* h.^2.*a.^2));
% 
% E_inf   = (1-i).*A1.*delta .* exp(-i.*k.*r) .* ...
%     (-(1./k.*r.^2) + 3*i./(k.^2.*r.^3) + 3./(k.^3.*r.^4));
% 
% 
% stimFlw = i.*omega.*c.dispAmp;
% bunFlw  = stimFlw .* (1-exp(-c.bunHeight.*(1+i)./delta));
% 
% S_bun       = c.M(1,:) ./ bunFlw;
% S_stimFlw   = bunFlw ./ stimFlw;
% S_sphere    = E_inf; 
% 
% 
% 
% 
% figure
% subplot(2,1,1)
% loglog(f,abs(S_bun),'r-'); hold on
% loglog(f,abs(S_bun.*S_stimFlw),'g'); 
% loglog(f,abs(S_bun.*S_stimFlw.*E_inf),'b'); 
% modPlot
% ylabel('Sensitivity (s)')
% xlabel('Freq. (Hz)');
% subplot(2,1,2)
% semilogx(f,angle(S_bun)*180/pi,'r-'); hold on
% semilogx(f,angle(S_bun.*S_stimFlw)*180/pi,'g');
% semilogx(f,angle(S_bun.*S_stimFlw.*E_inf)*180/pi,'b'); 
% modPlot
% ylabel('Phase (deg)');
% xlabel('Freq. (Hz)');







figure(f1)
subplot(2,1,1)
ylabel('Sensitivity (s)')
xlabel('Freq. (Hz)');
modPlot
subplot(2,1,2)
ylabel('Phase (deg)');
xlabel('Freq. (Hz)');
modPlot
    


function modPlot
set(gca,'TickDir','out')
set(gca,'TickLength',[.02 .02])

      
function c = c_default_theo
%Parameters for all anlayses
c.freqs         = [10.^linspace(-1,2,100)]';
c.numHeights    = 1;
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
c.numHairs      = 10;


function c = speedofSound(T)
%Gives speed of sound in m/s for temperature in degrees C
% From Bilaniuk & Wong (1993 & 1996)
c =   1.40238744 * 10^3 + 5.03836171 * T - 5.81172916 * 10^-2 * T^2 ...
    + 3.34638117 * 10^-4 * T^3 - 1.48259672 * 10^-6 * T^4 ...
    + 3.16585020 * 10^-9 * T^5;