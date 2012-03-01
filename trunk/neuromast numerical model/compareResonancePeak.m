function compareResonancePeak
% Runs the numerical model of the neuromast


c           = c_default_theo;
h           = c.bunHeight;
relativeTo  = 'freestream spd';
%relativeTo  = 'hair bundle spd';

a1          = mean([c.baseDiameter c.midDiameter])/2;
I_base      = (pi / 64) * (2 * a1)^4;
EI_mid      = (I_base * c.E_matrix) + (c.EI_kino * c.numHairs);


c.EIspecial = EI_mid / 1;
c1          = numerical_twopart(c,'torsion spring','special');
c1          = calcFreqResp(c1,h,relativeTo);


%c.EIspecial = EI_mid / 5;
%c2          = numerical_twopart(c,'torsion spring','one part');
%c2          = calcFreqResp(c2,h,relativeTo);

c.EIspecial = EI_mid / 5;
c2          = numerical_twopart(c,'torsion spring','special');
c2          = calcFreqResp(c2,h,relativeTo);


c.EIspecial = EI_mid / 10;
c3          = numerical_twopart(c,'torsion spring','special');
c3          = calcFreqResp(c3,h,relativeTo);


c.EIspecial = EI_mid / 15;
c4          = numerical_twopart(c,'torsion spring','special');
c4          = calcFreqResp(c4,h,relativeTo);









figure;
clrs = [1 0 0; .75 0 .75;0 0 0;0 .75 .75];
lins = {{'-'},{'-'},{'-'},{'-'}};
plotFreqResp({c1,c2,c3,c4}, ...
    clrs, lins,'1','2','3','4');
subplot(2,1,1)
title(['Bode for height of ' num2str(h) ' m']);



return



function c = c_default_theo
%Parameters for all anlayses
c.freqs         = [10.^linspace(-1,3,100)]';
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


function U = bLayer(x,omega,mu,rho,D0)

delta   = (2*mu./(rho.*omega)).^0.5;
disp    = D0 .* (1 - exp(-x.*(1+i) ./ delta));
U       = disp.*i.*omega;





