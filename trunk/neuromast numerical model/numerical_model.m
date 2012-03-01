function numerical_model
% Runs the numerical model of the neuromast

% 
% %Data from morphometric measurements (replace with real data)
% c.baseDiameter 	= 8.8e-6;
% c.midDiameter 	= 7.2e-6;
% c.kinoHeight 	= 29.7e-6;
% c.cupHeight     = 45e-6;
% c.numHairs      = 7;
% 
% %Parameters for the anlaysis
% %c.freqs         = 119.5792.* units.T; %[10.^linspace(2.077,2.078,300)]' .* units.T;
% c.freqs         = [10.^linspace(-3,3,100)]';
% c.numHeights    = 100;
% c.bunHeight     = 5.3e-6; %From Dinklo, 2005
% c.dispAmp       = .3 * 10^-6; %m
% c.E_matrix      = 31; %31 Pa
% c.EI_kino       = 2e-21; % 2e-21 N m^2
% c.bundleStiff   = 2.9e-14; %Nm/rad (van Netten, 1989) 
% c.rho           = 998; %998 kg m^-3
% c.mu            = 1.002e-3; %1.002e-3 Pa s

c   = c_default;
%c.dispAmp       = 100 * 10^-6; %m
%c.bundleStiff   = 10^5 * 2.925e-14;
%c.E_matrix      = 100;
%c.kinoHeight 	= c.cupHeight-.001e-6;
%c.midDiameter 	= 8.88e-6;
%c.freqs =5;
%c.numHairs      = 1;

%c.EI_kino       =  10 * 2e-21;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = numerical_twopart(c,'torsion spring','two part');


%c.dispAmp = 10^5*c.dispAmp 

% anlyze model
%h = c.bunHeight;
h = 5.3e-6;
c = calcFreqResp(c,h);




%compareToSimple(c)





%figure;
clr = [1 .2 1];
clr = [1 0 0];
plotFreqResp({c},clr)

 %figure;
 %visDeflect(c,1)
 %figure;
 %visDeflect(c,10)
% figure;
% visDeflect(c,100)
% figure;

%visForces(c,[.1 10^-.5 100])

return

 disp(' ')
 disp(['This neuromast has a sensitivity of ' num2str(c.peak_sense) ' s'])
 disp(['This peak ocurs at ' num2str(c.peak_freq) ' Hz'])
 disp(' ')

 
 
 
function c = c_default
%Parameters for all anlayses
 c.freqs         = [10.^linspace(-2,3,200)]';
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
c.kinoHeight 	= 15e-6;
c.cupHeight     = 45e-6;
c.numHairs      = 11;

