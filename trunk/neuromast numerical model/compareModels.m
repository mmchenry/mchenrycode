function compareModels(action)
% Runs the numerical model of the neuromast

if nargin < 1
    action = 'all models';
end

switch action
    case 'all models'
        c = c_default_theo;
        %h = c.cupHeight;
        h = c.bunHeight;

        %MODEL 1: TWO-PART, FIXED BASE, NUMERICAL
        %-------------------------------------------------------------------------
        %The model has realistic parameter values


        c   = c_default_theo;
        c1  = numerical_twopart(c,'fixed');
        c1  = calcFreqResp(c1,h,'freestream spd');

        clear c
        

        %MODEL 2: TWO-PART, TORSION SPRING, NUMERICAL
        %-------------------------------------------------------------------------
        %The model has realistic parameter values


        c   = c_default_theo;
        %c.dispAmp       = 5 * 10^-4;
        c2  = numerical_twopart(c,'torsion spring');
        %c2  = calcFreqResp(c2,h);
        c2  = calcFreqResp(c2,h,'freestream spd');
        %visForces(c2,[.1 10 100 1000])
        clear c


        %MODEL 3: TWO-PART, LINEAR SPRING, NUMERICAL
        %-------------------------------------------------------------------------
        %The model has realistic parameter values

        c   = c_default_theo;
        %c.linStiff      = 10^-1 .* 0.13 * 10^-3;
        c3  = numerical_twopart(c,'linear spring');
        c3  = calcFreqResp(c3,h,'freestream spd');

        clear c


        %MODEL 4: ONE-PART, TORSION SPRING, LOW FREQ ANALYTICAL
        %-------------------------------------------------------------------------

        %Data from morphometric measurements (based on stiffness paper)
        c       = c_default_theo;
        %c.freqs = [10.^linspace(-5,-2,20)]';

        c4      = runAnalytical(c,'low');
        c4      = calcFreqResp(c4,h,'freestream spd');

        clear c


        %MODEL 5: ONE-PART, TORSION SPRING, NUMERICAL
        %-------------------------------------------------------------------------
        % Makes the tip of the second part very small, extends first part

        %Data from morphometric measurements (based on stiffness paper)
        c   = c_default_theo;
        c.baseDiameter 	= 8.88e-6;
        c.midDiameter 	= 8.88e-6;
        %c.kinoHeight 	= 44.999e-6;
        %c.kinoHeight 	= .00001e-6;
        c5 = numerical_twopart(c,'torsion spring','one part');
        c5 = calcFreqResp(c5,h,'freestream spd');

        clear c


        %MODEL 6: ONE-PART, TORSION SPRING, HIGH FREQ ANALYTICAL
        %-------------------------------------------------------------------------

        %Data from morphometric measurements (based on stiffness paper)
        c       = c_default_theo;
        %c.freqs = [10.^linspace(-5,-2,20)]';

        c6      = runAnalytical(c,'high');
        c6      = calcFreqResp(c6,h,'freestream spd');

        clear c



        % figure;
        % clrs = [1 0 0];
        % lins = {{'-'}};
        % plotFreqResp({c5},...
        %             clrs, lins)
        % subplot(2,1,1)
        % title(['Bode for height of ' num2str(h) ' m']);


        figure;
        clrs = [1 0 0; .75 0 .75;0 0 0;0 .75 .75;0 0 1];
        lins = {{'-'},{'-'},{'--'},{'-'},{'-'}};
        plotFreqResp({c1,c2,c3,c4,c5}, ...
            clrs, lins,'2 pins','1 pin','wheels','low f','1 part');
        subplot(2,1,1)
        title(['Bode for height of ' num2str(h) ' m']);

        
        
    case 'different params'  
        
        c   = c_default_theo;
        h = c.bunHeight;
        c1  = numerical_twopart(c,'torsion spring');
        c1  = calcFreqResp(c1,h,'freestream spd');
        clear c
        
        
        
         c   = c_default_theo;
         hDensity = 0.13; %hair cells / micron^2
 
        c.baseDiameter 	= 200e-6;
        rat             = c.midDiameter/c.baseDiameter;
        c.midDiameter   = rat * c.baseDiameter;
        c.numHairs      = hDensity * (pi*((c.baseDiameter*10^6)/2)^2); 
        
        c2  = numerical_twopart(c,'fixed');
        c2  = calcFreqResp(c2,h,'freestream spd');

        clear c
        
        
        
        c               = c_default_theo;
        AR              = c.cupHeight/c.baseDiameter;
        c.baseDiameter 	= 50e-6;
        rat             = c.midDiameter/c.baseDiameter;
        c.midDiameter   = rat * c.baseDiameter;
        c.numHairs      = hDensity * (pi*((c.baseDiameter*10^6)/2)^2);
        c.cupHeight     = AR .* c.baseDiameter;
        
        c3  = numerical_twopart(c,'fixed');
        c3  = calcFreqResp(c3,h,'freestream spd');

        clear c
        

        
        figure;
        clrs = [1 0 0; .75 0 .75;0 .75 .75];
        lins = {{'-'},{'-'},{'-'}};
        plotFreqResp({c1,c2,c3}, clrs, lins,' ',' ',' ');
        subplot(2,1,1)
        title(['Bode for height of ' num2str(h) ' m']);
        
        
        
    case 'structure and bLayer'
        c = c_default_theo;
        h = c.bunHeight;
        h = 9e-6;
        
        %c.numHairs      = 1;
        %c.baseDiameter 	= 10;
        
        %MODEL 1: STRUCTURE [TWO-PART, TORSION SPRING, NUMERICAL]
        %------------------------------------------------------------------
        c1  = numerical_twopart(c,'torsion spring');
        c1  = calcFreqResp(c1,h,'hair bundle spd');
        
        %MODEL 2: FLUID 
        %------------------------------------------------------------------
        c2  = runAnalytical(c,'bLayer');
        c2  = calcFreqResp(c2,h,'freestream spd');

        %MODEL 3: RELATIVE TO FREESTREAM
        %------------------------------------------------------------------
        c3  = numerical_twopart(c,'torsion spring');
        c3  = calcFreqResp(c3,h,'freestream spd');
        
        
%         figure;
%         subplot(3,1,1)
%         loglog(c1.freqs,c1.sensitivity,'k-')
%         subplot(3,1,2)
%         loglog(c2.freqs,c2.sensitivity,'k-')
%         subplot(3,1,3)
%         loglog(c3.freqs,c3.sensitivity,'k-')
        
        figure;
        subplot(2,1,1)
        loglog(c1.freqs,c1.sensitivity,'k-');hold on
        set(gca,'TickDir','out')
        set(gca,'TickLength',[.02 .02])
        title('relative to bundle');
        subplot(2,1,2)
        semilogx(c1.freqs,c1.phase,'k-');hold on
        
        subplot(2,1,1)
        %plotyy(c3.freqs,c3.sensitivity,c2.freqs,c2.sensitivity,'loglog')
        loglog(c3.freqs,c3.sensitivity)
        set(gca,'TickDir','out')
        set(gca,'TickLength',[.02 .02])
        title('relative to freestream');
        subplot(2,1,2)
        semilogx(c3.freqs,c3.phase);
        
        figure;
        subplot(2,1,1)
        loglog(c2.freqs,c2.sensitivity,'r');
        set(gca,'TickDir','out')
        set(gca,'TickLength',[.02 .02])
        title('fluid');
        subplot(2,1,2)
        semilogx(c2.freqs,c2.phase);
%         
%         figure;
%         clrs = [1 0 0; .75 0 .75];
%         lins = {{'-'},{'-'}};
%         plotFreqResp({c1,c3},clrs, lins,' ',' ');
%         subplot(2,1,1)
%         title([' ']);
%         
%         figure;
%         clrs = [1 0 0; .75 0 .75;0 0 0];
%         lins = {{'-'},{'-'},{'-'}};
%         plotFreqResp({c1,c2,c3},clrs, lins,' ',' ',' ');
%         subplot(2,1,1)
%         title([' ']);

    case 'just bLayer'
        c = c_default_theo;
        x = linspace(1e-8,c.cupHeight,100);
        c.freqs = [10.^linspace(-1,5,200)]';
        %U = bLayer(x,omega,c.mu,c.rho,D0);
        
        
        %BODE PLOTS
        %------------------------------------------------------------------
        figure;
        
        omega   = 2*pi.*c.freqs;
        delta   = (2*c.mu./(c.rho.*omega)).^0.5;
        

        c2  = runAnalytical(c,'bLayer');
        xs = [5:10:35] .* 10^-6;
        for j=1:length(xs)
            relAmp  = (1 - exp(-xs(j).*(1+i) ./ delta));
            subplot(2,1,1)
            loglog(c.freqs,abs(relAmp),'r');
            
            hold on
            set(gca,'TickDir','out')
            set(gca,'TickLength',[.02 .02])
            title('fluid');
            subplot(2,1,2);
            semilogx(c.freqs,angle(relAmp)*180/pi);
            set(gca,'TickDir','out')
            set(gca,'TickLength',[.02 .02])
            hold on
        end
        
        clear omega delta xs relAmp
        
        
        
        %VS HEIGHT
        %------------------------------------------------------------------
        figure;

        L1  = c.kinoHeight;
        freqs = 10.^[0:4];
        
        for j=1:length(freqs)
            omega   = 2*pi.*freqs(j);
            delta   = (2*c.mu./(c.rho.*omega)).^0.5;
            relAmp  = (1 - exp(-x.*(1+i) ./ delta));
            
            clr = 'b';
            
            subplot(1,2,1)
            h = plot(abs(relAmp),10^6 .* x,'-');hold on
            temp = get(gca,'XTick');
            h2 = plot([min(temp) max(temp)],10^6 .*[L1 L1],'r-');
            ylabel('Height (micrn)');
            xlabel('Amp (micrn)')
            axis square
            hold on
            set(gca,'TickDir','out')
            set(gca,'TickLength',[.02 .02])
            
            subplot(1,2,2)
            h = plot(angle(relAmp)*180/pi,10^6 .* x,'-');hold on
            temp = get(gca,'XTick');
            h2 = plot([min(temp) max(temp)],10^6 .*[L1 L1],'r-');
            xlabel('Phase (deg)');
            ylabel('Height (micrn)');
            set(gca,'TickDir','out')
            set(gca,'TickLength',[.02 .02])
            axis square
            hold on
        end
        subplot(1,2,1)
        tmp=get(gca,'XLim');
        set(gca,'XLim',[tmp(1)-range(tmp)/20 tmp(2)]); clear tmp
        

end


% figure;
% visDeflect(c,1)
% figure;
% visDeflect(c,10)
% figure;
% visDeflect(c,100)
% figure;
%visForces(c,[.01 10 100 1000])






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


function c = c_default_will
%Mean morphometric values from Will's data
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
c.baseDiameter 	= 9.8e-6;
c.midDiameter 	= 10.9e-6 ;
c.kinoHeight 	= 21.7e-6;
c.cupHeight     = 35.7e-6;
c.numHairs      = 9; %Harris_2003


function c = runAnalytical(c,mType)

%Define general parameters
E_matrix    = c.E_matrix;
EI_kino     = c.EI_kino;
bundleStiff = c.bundleStiff;
rho         = c.rho;
mu          = c.mu; 
D0          = c.dispAmp;
L           = c.cupHeight;
a           = mean([c.baseDiameter c.midDiameter])/2;
freqs       = c.freqs;

I           = (pi / 64) * (2 * a)^4;
EI          = (I * E_matrix) + (EI_kino * c.numHairs);
baseStiff   = bundleStiff .* c.numHairs;
linStiff    = c.linStiff .* c.numHairs;

num1        = round(c.numHeights .* c.kinoHeight / c.cupHeight);
x1          = [linspace(0,c.kinoHeight, num1)]';
x1          = x1(2:end);
x2          = [linspace(c.kinoHeight,c.cupHeight,c.numHeights-num1+1)]';
x2          = x2 - x2(1);
c.heights   = [x1(1:end-1);x2+x1(end)];

for ii = 1:length(freqs)
    omega       = 2*pi*freqs(ii);
    if strcmp(mType,'low')
        c.M(:,ii)   = simpleBeam_low_2pin(c.heights,omega,mu,rho,L,D0...
                                        ,a,EI,baseStiff);
    elseif strcmp(mType,'high') 
        c.M(:,ii)   = simpleBeam_high(c.heights,omega,mu,rho,L,D0...
                                        ,a,EI,linStiff);
    elseif strcmp(mType,'bLayer') 
        c.M(:,ii)   = bLayer(c.heights,omega,mu,rho,D0);
    else
        error('incorrect model type');
    end
end

function v = simpleBeam_low(x,omega,mu,rho,L,D0,a,EI,kn)
%one-part model for low freq - all fluid forces

delta   = (2 * mu / rho / omega)^0.5;
s       = a ./ 2^.5 / delta;
g       = 0.577 + log(s);
G       = -g / (g^2 + (pi/4)^2);
d       = D0 * (1 + i) * x ./ delta;
b       = omega * pi * (2*rho*a^2*omega - 4*mu*G*i - pi*mu*G/g);
%b       = omega * pi * (- 4*mu*G*i);
v       = -(1/120 + i/120) * b*D0.*x .* ...
            (40*EI*L^3 + kn.*x .* (20*L^3 - 10*L^2.*x + x.^3)) ./ ...
            (EI * kn *delta);

function v = simpleBeam_low_2pin(x,omega,mu,rho,L,D0,a,EI,kn)
%one-part model for low freq - all fluid forces - anchored base

delta   = (2 * mu / rho / omega)^0.5;
s       = a ./ 2^.5 / delta;
g       = 0.577 + log(s);
G       = -g / (g^2 + (pi/4)^2);
d       = D0 * (1 + i) * x ./ delta;
b       = omega * pi * (2*rho*a^2*omega - 4*mu*G*i - pi*mu*G/g);
v       = -(1/120 + i/120) * b*D0.*x.^2 .* ...
            (20*L^3-10*L^2.*x+x.^3) ./ (EI*delta);        
        
% function v = simpleBeam_high(x,omega,mu,rho,L,D0,a,EI,kn)
% %one-part model for low freq - all fluid forces
% 
% delta   = (2 * mu / rho / omega)^0.5;
% s       = a ./ 2^.5 / delta;
% g       = 0.577 + log(s);
% G       = -g / (g^2 + (pi/4)^2);
% d       = D0 * (1 + i) * x ./ delta;
% b       = omega * pi * (2*rho*a^2*omega - 4*mu*G*i - pi*mu*G/g);
% %b       = omega * pi * (2*rho*a^2*omega - pi*mu*G/g);
% v       = - b*D0.*x .* (12*EI*L^2 + kn.*x .* (6*L^2 - 4*L.*x + x.^2)) ./ ...
%             (24.*EI.*kn);

function v = simpleBeam_high(x,omega,mu,rho,L,D0,a,EI,k)
%one-part model for low freq - all fluid forces

delta   = (2 * mu / rho / omega)^0.5;
s       = a ./ 2^.5 / delta;
g       = 0.577 + log(s);
G       = -g / (g^2 + (pi/4)^2);
d       = D0 * (1 + i) * x ./ delta;
%b       = omega * pi * (2*rho*a^2*omega - 4*mu*G*i - pi*mu*G/g);
b       = omega * pi * (2*rho*a^2*omega - pi*mu*G/g);
v       = - b*D0.*L .*(6.*EI + k.*x.^2 .* (-3.*L + x)) ./ ...
            (6*EI*k);

function U = bLayer(x,omega,mu,rho,D0)

delta   = (2*mu./(rho.*omega)).^0.5;
disp    = D0 .* (1 - exp(-x.*(1+i) ./ delta));
U       = disp.*i.*omega;









