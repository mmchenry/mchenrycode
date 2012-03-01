function numerical_dinklo
% Runs the numerical model of the neuromast

load dinkloData



h           = 29e-6;

%h=5e-6

phOffset    = 0;
c           = c_default;
%h           = c.bunHeight
%c.E_matrix  = 10; 
%c.EI_kino   = 2.0e-21;
c.cupHeight = 32e-6;
%c.kinoHeight= 18e-6;
% c.baseDiameter 	= 8.88e-6;
% c.midDiameter 	= 6e-6 ;
%c.numHairs  = 11;

% h           = 5e-6;
% phOffset    = 65;
% c           = c_default;
% c.E_matrix  = 20; 
% c.EI_kino   = 2e-21;
% c.cupHeight = 40e-6;
% c.numHairs  = 7;


c   = numerical_twopart(c,'torsion spring');
c   = calcFreqResp(c,h,'freestream spd');



figure;

%subplot(2,1,1)

h = loglog((c.freqs),(c.sensitivity).*10.^0,'-k');
hold on
loglog(dinklo.freq,dinklo.sense.*10.^-3,'or');
%grid on;
ylabel('Sensitivity (\mu{s})')
xlabel('Freq. (Hz)');
axis square
%set(h,'Color',clr(ii,:))
%set(h,'MarkerEdgeColor',clr(ii,:))
set(gca,'TickDir','out')
    set(gca,'TickLength',[.02 .02])
grid on
    return

subplot(2,1,2)

h = semilogx((c.freqs),c.phase,'-k');
hold on
semilogx(dinklo.freq,dinklo.phase + phOffset,'or');
ylabel('Phase (deg)');
xlabel('Freq. (Hz)');
%grid on
set(gca,'TickDir','out')
    set(gca,'TickLength',[.02 .02])
axis square
%set(h,'Color',clr(ii,:))


%figure;
%visDeflect(c,1)




% set(h,'MarkerSize',3)
% set(h,'MarkerEdgeColor',clr)
% set(h,'MarkerFaceColor',clr)
 
 
 
 
 
 
 
 
function c = c_default
%Parameters for all anlayses
 c.freqs         = [10.^linspace(0,3,200)]';
 c.numHeights    = 50;  
 c.bunHeight     = 5.3e-6; %From Dinklo, 2005
 c.dispAmp       = 10 * 10^-6; %m
 c.E_matrix      = 31; %31 Pa
 c.EI_kino       = 2.0e-21; % 2e-21 N m^2
 c.bundleStiff   = 2.925e-14; %Nm/rad (van Netten & Kroese, 1987) 
 c.linStiff      = 0.13 * 10^-3; %N/m (van Netten & Kroese, 1987) 
 c.rho           = 998; %998 kg m^-3
 c.mu            = 1.002e-3; %1.002e-3 Pa s

%Data from morphometric measurements (based on stiffness paper)
c.baseDiameter 	= 8.88e-6;
c.midDiameter 	= 7.2e-6 ;
c.kinoHeight 	= 16e-6;
c.cupHeight     = 45e-6;
c.numHairs      = 11;