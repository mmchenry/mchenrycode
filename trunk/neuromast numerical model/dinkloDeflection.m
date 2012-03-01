function dinkloDeflection


%Comparison with Dinklo's experiments (vs. height at 44 Hz)
%-------------------------------------------------------------------------




%Model parameters to match Theo's experiments


%c.kinoHeight 	= 29.7e-6;

c           = c_default;
%c.kinoHeight= 18e-6;
c.dispAmp       = .26 * 10^-6; %m


c.freqs         = 44;
%c.numHairs  = 11;

c.cupHeight     = 38e-6;

%c.numHeights    = 1;

c  = numerical_twopart(c,'torsion spring');

vVal    = c.M(:,1);
h       = c.heights;

%Boundary layer flow
omega   = 2 * pi .* c.freqs;
delta   = ( (2 .* c.mu) ./ (c.rho .* omega) ).^.5;
freeStrm= c.dispAmp .* i .* omega;

d       =  freeStrm .* (1 - exp( -(1+i) .* h ./ delta)); 
hFLow   = freeStrm .* (1 - exp( -(1+i) .* c.bunHeight ./ delta)); 


% Get Theo's data:

load('theoData')

time        = theo.time.all.*1000;
sRate       = 1 / mean( diff( time ) );
p1Y         = theo.pt1.y;
p2Y         = theo.pt2.y;
p3Y         = theo.pt3.y;

p1FFT       = fft( theo.pt1.xAll );
p2FFT       = fft( theo.pt2.xAll );
p3FFT       = fft( theo.pt3.xAll );

p1AMP       = abs( 2.* p1FFT ./ length(p1FFT) );
p2AMP       = abs( 2.* p2FFT ./ length(p2FFT) );
p3AMP       = abs( 2.* p3FFT ./ length(p3FFT) );

p1PH        = angle( p1FFT )  .* (180 / pi);
p2PH        = angle( p2FFT  )  .* (180 / pi);
p3PH        = angle( p3FFT  )  .* (180 / pi);

f           = (0:length(p3FFT)-1)'*sRate/length(p3FFT);
iFreq       = 1:round(length(f)/2);
iPeak       = min( find( p1AMP(iFreq) == max( p1AMP(iFreq) )) );
p1AMP       = p1AMP( iPeak );
p2AMP       = p2AMP( iPeak );
p3AMP       = p3AMP( iPeak );
p1PH        = p1PH( iPeak );
p2PH        = p2PH( iPeak );
p3PH        = p3PH( iPeak );

nudge =  -42.9442 - p3PH;

figure;

yBoost = 0*2e-6;
subplot(1,2,1)

plot(abs(vVal), h, 'k-');
hold on

plot( [ p1AMP p2AMP p3AMP ] , [ p1Y p2Y p3Y ] + yBoost , 'or')
modPlot
%h3 = plot(abs(d),h);
%set(h3,'Color',.8.*[1 1 1])

ylabel('Height (m)')
xlabel('Amplitude (m)');
%axis equal


subplot(1,2,2)

plot(  180*(angle(vVal./freeStrm))/pi, h, 'k-');
hold on

plot( [ p1PH p2PH p3PH ] + nudge, [ p1Y p2Y p3Y ] + yBoost, 'or')
modPlot
%h3 = plot((180/pi)* (angle(d./freeStrm))- 45 , h);
%set(h3,'Color',.8.*[1 1 1])

xlabel('Phase (deg)');
%axis square
legend('cupula','theo data');



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
c.kinoHeight 	= 15e-6;
c.cupHeight     = 45e-6;
c.numHairs      = 11;


function modPlot
set(gca,'TickDir','out')
set(gca,'TickLength',[.02 .02])
%set(gca,'YLim',[0 max(heights)])
%set(gca,'YTick',[0:1e-5:max(heights)])
tmp=get(gca,'XLim');
set(gca,'XLim',[tmp(1)-range(tmp)/20 tmp(2)]); clear tmp
tmp=get(gca,'YLim');
set(gca,'YLim',[tmp(1)-range(tmp)/20 tmp(2)]); clear tmp
%axis square

