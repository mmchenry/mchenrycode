function d = sensitivityAnalysis_varyParam


wrt = 'freestream';

c           = c_default_theo;
numCols     = 4;
numRows     = 4;
senseLim    = [0 .0003];
freqLim     = [0 25];

% if strcmp(wrt,'freestream')
%     c.freqs         = [10.^linspace(-0.5,1.5,500)]';
% end


% % Run  'num hairs'
% 
% vals    = round(linspace(2,100,100));
% exVals  = [8 20];
% 
% cT      = c;
% 
% for j = 1:length(vals)
%     cT.numHairs = vals(j);
% 
%     cT  = numerical_twopart(cT,'torsion spring');
%     cT  = calcBunResp(cT,wrt);
%     [c_sense,c_cut] = calcPerformance(cT,wrt);
% 
%     dHairs.val(j)       = vals(j);
%     dHairs.Dsense(j)    = c_sense;
%     dHairs.Dcut(j)      = c_cut;
% 
% end
% 
% cE1             = c;
% cE1.numHairs    = exVals(1);
% cE1.freqs       = [10.^linspace(-2,3,500)]';
% cE1             = numerical_twopart(cE1,'torsion spring');
% cE1             = calcBunResp(cE1,wrt);
% 
% cE2             = c;
% cE2.numHairs    = exVals(2);
% cE2.freqs       = cE1.freqs;
% cE2             = numerical_twopart(cE2,'torsion spring');
% cE2             = calcBunResp(cE2,wrt);
% 
% 
% 
 f1 = figure;
% subplot(numRows,numCols,1)
% plot(dHairs.val,dHairs.Dsense);hold on
% plot([exVals(1) exVals(1)],[senseLim(1) senseLim(2)],'g--')
% plot([exVals(2) exVals(2)],[senseLim(1) senseLim(2)],'r--')
% set(gca,'YLim',senseLim)
% modPlot
% ylabel('max sensitivity (s)')
% xlabel('num hair cells')
% 
% subplot(numRows,numCols,numCols+1)
% plot(dHairs.val,dHairs.Dcut);hold on
% plot([exVals(1) exVals(1)],[freqLim(1) freqLim(2)],'g--')
% plot([exVals(2) exVals(2)],[freqLim(1) freqLim(2)],'r--')
% %set(gca,'YLim',freqLim)
% modPlot
% ylabel('cut freq(Hz)')
% xlabel('num hair cells')
% 
 f2=figure;
% subplot(2,numCols,1)
% loglog(cE1.freqs,cE1.sensitivity,'g',cE2.freqs,cE2.sensitivity,'r--');
% set(gca,'TickDir','out')
% set(gca,'TickLength',[.04 .04])
% set(gca,'XLim',[0.01 1000])
% set(gca,'XTick',[10.^[-1:1:3]])
% set(gca,'YLim',10.^[-6 -3])
% set(gca,'YTick',[10.^[-6:1:-3]])
% 
% clear cT cE1 cE2 dDai
% 
% 
% 
% % Run  'kino height'
% vals    = linspace(5,44.9,50) .* 10^-6;
% exVals  = [16 30].*10^-6;
% cT      = c;
% 
% for j = 1:length(vals)
%     cT.kinoHeight = vals(j);
% 
%     cT  = numerical_twopart(cT,'torsion spring');
%     cT  = calcBunResp(cT,wrt);
%     [c_sense,c_cut] = calcPerformance(cT,wrt);
% 
%     dKino.val(j)       = vals(j);
%     dKino.Dsense(j)    = c_sense;
%     dKino.Dcut(j)      = c_cut;
% 
% end
% 
% 
% cE1             = c;
% cE1.kinoHeight  = exVals(1);
% cE1.freqs       = [10.^linspace(-2,3,500)]';
% cE1             = numerical_twopart(cE1,'torsion spring');
% cE1             = calcBunResp(cE1,wrt);
% 
% cE2             = c;
% cE2.kinoHeight  = exVals(2);
% cE2.freqs       = cE1.freqs;
% cE2             = numerical_twopart(cE2,'torsion spring');
% cE2             = calcBunResp(cE2,wrt);
% 
% 
% 
% figure(f1)
% subplot(numRows,numCols,2)
% plot(dKino.val.*10^6,dKino.Dsense);hold on
% plot([exVals(1) exVals(1)].*10^6,[senseLim(1) senseLim(2)],'g--')
% plot([exVals(2) exVals(2)].*10^6,[senseLim(1) senseLim(2)],'r--')
% set(gca,'YLim',senseLim)
% modPlot
% ylabel('max sensitivity (s)')
% xlabel('kino height')
% subplot(numRows,numCols,numCols+2)
% plot(dKino.val.*10^6,dKino.Dcut);hold on
% plot([exVals(1) exVals(1)].*10^6,[freqLim(1) freqLim(2)],'g--')
% plot([exVals(2) exVals(2)].*10^6,[freqLim(1) freqLim(2)],'r--')
% set(gca,'YLim',freqLim)
% modPlot
% ylabel('cut freq(Hz)')
% xlabel('kino height')
% 
% figure(f2);
% subplot(2,numCols,2)
% loglog(cE1.freqs,cE1.sensitivity,'g',cE2.freqs,cE2.sensitivity,'r--');
% set(gca,'TickDir','out')
% set(gca,'TickLength',[.04 .04])
% set(gca,'XLim',[0.01 1000])
% set(gca,'XTick',[10.^[-1:1:3]])
% set(gca,'YLim',10.^[-6 -3])
% set(gca,'YTick',[10.^[-6:1:-3]])
% 
% clear cT cE1 cE2 dDai


% % Run  'cup height'
% vals    = linspace(16,100,70) .* 10^-6;
% exVals  = [30 90].*10^-6;
% cT      = c;
% cT.kinoHeight = 15e-6;
% 
% for j = 1:length(vals)
%     cT.cupHeight = vals(j);
% 
%     cT  = numerical_twopart(cT,'torsion spring');
%     cT  = calcBunResp(cT,wrt);
%     [c_sense,c_cut] = calcPerformance(cT,wrt);
% 
%     dCup.val(j)       = vals(j);
%     dCup.Dsense(j)    = c_sense;
%     dCup.Dcut(j)      = c_cut;
% 
% end
% 
% 
% cE1             = c;
% cE1.cupHeight   = exVals(1);
% cE1.freqs       = [10.^linspace(-2,3,500)]';
% cE1             = numerical_twopart(cE1,'torsion spring');
% cE1             = calcBunResp(cE1,wrt);
% 
% cE2             = c;
% cE2.cupHeight   = exVals(2);
% cE2.freqs       = cE1.freqs;
% cE2             = numerical_twopart(cE2,'torsion spring');
% cE2             = calcBunResp(cE2,wrt);
% 
% 
% 
% figure(f1)
% subplot(numRows,numCols,3)
% plot(dCup.val.*10^6,dCup.Dsense); hold on
% plot([exVals(1) exVals(1)].*10^6,[senseLim(1) senseLim(2)],'g--')
% plot([exVals(2) exVals(2)].*10^6,[senseLim(1) senseLim(2)],'r--')
% set(gca,'YLim',senseLim)
% modPlot
% ylabel('max sensitivity (s)')
% xlabel('cupula height')
% subplot(numRows,numCols,numCols+3)
% plot(dCup.val.*10^6,dCup.Dcut);hold on
% plot([exVals(1) exVals(1)].*10^6,[freqLim(1) freqLim(2)],'g--')
% plot([exVals(2) exVals(2)].*10^6,[freqLim(1) freqLim(2)],'r--')
% set(gca,'YLim',freqLim)
% modPlot
% ylabel('cut freq(Hz)')
% xlabel('cupula height')
% 
% figure(f2);
% subplot(2,numCols,3)
% loglog(cE1.freqs,cE1.sensitivity,'g',cE2.freqs,cE2.sensitivity,'r--');
% set(gca,'TickDir','out')
% set(gca,'TickLength',[.04 .04])
% set(gca,'XLim',[0.01 1000])
% set(gca,'XTick',[10.^[-1:1:3]])
% set(gca,'YLim',10.^[-6 -3])
% set(gca,'YTick',[10.^[-6:1:-3]])
% 
% clear cT cE1 cE2 dDai
% 
% 
% 
% % Run  'cup diameter'
% vals    = linspace(8,30,100) .* 10^-6;
% exVals  = [10 25].*10^-6;
% cT      = c;
% rat = cT.midDiameter/cT.baseDiameter;
% 
% for j = 1:length(vals)
% 
%     cT.baseDiameter = vals(j);
%     cT.midDiameter = rat * vals(j);
% 
%     cT  = numerical_twopart(cT,'torsion spring');
%     cT  = calcBunResp(cT,wrt);
%     [c_sense,c_cut] = calcPerformance(cT,wrt);
% 
%     dDai.val(j)       = vals(j);
%     dDai.Dsense(j)    = c_sense;
%     dDai.Dcut(j)      = c_cut;
% 
% end
% 
% cE1             = c;
% cE1.baseDiameter= exVals(1);
% cE1.midDiameter = rat * exVals(1);
% cE1.freqs       = [10.^linspace(-2,3,500)]';
% cE1             = numerical_twopart(cE1,'torsion spring');
% cE1             = calcBunResp(cE1,wrt);
% 
% cE2             = c;
% cE2.baseDiameter= exVals(2);
% cE2.midDiameter = rat * exVals(2);
% cE2.freqs       = cE1.freqs;
% cE2             = numerical_twopart(cE2,'torsion spring');
% cE2             = calcBunResp(cE2,wrt);
% 
% figure(f1)
% subplot(numRows,numCols,4)
% plot(dDai.val.*10^6,dDai.Dsense); hold on
% plot([exVals(1) exVals(1)].*10^6,[senseLim(1) senseLim(2)],'g--')
% plot([exVals(2) exVals(2)].*10^6,[senseLim(1) senseLim(2)],'r--')
% set(gca,'YLim',senseLim)
% modPlot
% ylabel('max sensitivity (s)')
% xlabel('cupula diameter')
% subplot(numRows,numCols,numCols+4)
% plot(dDai.val.*10^6,dDai.Dcut);hold on
% plot([exVals(1) exVals(1)].*10^6,[freqLim(1) freqLim(2)],'g--')
% plot([exVals(2) exVals(2)].*10^6,[freqLim(1) freqLim(2)],'r--')
% set(gca,'YLim',freqLim)
% modPlot
% ylabel('cut freq(Hz)')
% xlabel('cupula diameter')
% 
% figure(f2);
% subplot(2,numCols,4)
% loglog(cE1.freqs,cE1.sensitivity,'g',cE2.freqs,cE2.sensitivity,'r--');
% set(gca,'TickDir','out')
% set(gca,'TickLength',[.04 .04])
% set(gca,'XLim',[0.01 1000])
% set(gca,'XTick',[10.^[-1:1:3]])
% set(gca,'YLim',10.^[-6 -3])
% set(gca,'YTick',[10.^[-6:1:-3]])
% 
% clear cT cE1 cE2 dDai


% Run  'E_matrix'
vals    = linspace(10,150,100);
exVals  = [20 110];
cT      = c;

for j = 1:length(vals)

    cT.E_matrix = vals(j);

    cT  = numerical_twopart(cT,'torsion spring');
    cT  = calcBunResp(cT,wrt);
    [c_sense,c_cut] = calcPerformance(cT,wrt);

    dDai.val(j)       = vals(j);
    dDai.Dsense(j)    = c_sense;
    dDai.Dcut(j)      = c_cut;

end

cE1             = c;
cE1.E_matrix   = exVals(1);
cE1.freqs       = [10.^linspace(-2,3,500)]';
cE1             = numerical_twopart(cE1,'torsion spring');
cE1             = calcBunResp(cE1,wrt);

cE2             = c;
cE2.E_matrix    = exVals(2);
cE2.freqs       = cE1.freqs;
cE2             = numerical_twopart(cE2,'torsion spring');
cE2             = calcBunResp(cE2,wrt);


figure(f1)
subplot(numRows,numCols,9)
plot(dDai.val,dDai.Dsense); hold on
plot([exVals(1) exVals(1)],[senseLim(1) senseLim(2)],'g--')
plot([exVals(2) exVals(2)],[senseLim(1) senseLim(2)],'r--')
set(gca,'YLim',senseLim)
modPlot
ylabel('max sensitivity (s)')
xlabel('E matrix')
subplot(numRows,numCols,numCols+9)
plot(dDai.val,dDai.Dcut);hold on
plot([exVals(1) exVals(1)],[freqLim(1) freqLim(2)],'g--')
plot([exVals(2) exVals(2)],[freqLim(1) freqLim(2)],'r--')
set(gca,'YLim',freqLim)
modPlot
ylabel('cut freq(Hz)')
xlabel('E matrix')


figure(f2);
subplot(2,numCols,5)
loglog(cE1.freqs,cE1.sensitivity,'g',cE2.freqs,cE2.sensitivity,'r--');
set(gca,'TickDir','out')
set(gca,'TickLength',[.04 .04])
set(gca,'XLim',[0.01 1000])
set(gca,'XTick',[10.^[-1:1:3]])
set(gca,'YLim',10.^[-6 -3])
set(gca,'YTick',[10.^[-6:1:-3]])

clear cT cE1 cE2 dDai



% Run  'bundle stiff'
vals    = linspace(.01e-14,5e-14,50);
exVals  = [1e-14 4e-14];
cT      = c;

for j = 1:length(vals)

    cT.bundleStiff = vals(j);

    cT  = numerical_twopart(cT,'torsion spring');
    cT  = calcBunResp(cT,wrt);
    [c_sense,c_cut] = calcPerformance(cT,wrt);

    dDai.val(j)       = vals(j);
    dDai.Dsense(j)    = c_sense;
    dDai.Dcut(j)      = c_cut;

end

cE1             = c;
cE1.bundleStiff  = exVals(1);
cE1.freqs       = [10.^linspace(-2,3,500)]';
cE1             = numerical_twopart(cE1,'torsion spring');
cE1             = calcBunResp(cE1,wrt);

cE2             = c;
cE2.bundleStiff  = exVals(2);
cE2.freqs       = cE1.freqs;
cE2             = numerical_twopart(cE2,'torsion spring');
cE2             = calcBunResp(cE2,wrt);


figure(f1)
subplot(numRows,numCols,10)
plot(dDai.val,dDai.Dsense); hold on
plot([exVals(1) exVals(1)],[senseLim(1) senseLim(2)],'g--')
plot([exVals(2) exVals(2)],[senseLim(1) senseLim(2)],'r--')
set(gca,'YLim',senseLim)
modPlot
ylabel('max sensitivity (s)')
xlabel('bundle stiff')
subplot(numRows,numCols,numCols+10)
plot(dDai.val,dDai.Dcut);hold on
plot([exVals(1) exVals(1)],[freqLim(1) freqLim(2)],'g--')
plot([exVals(2) exVals(2)],[freqLim(1) freqLim(2)],'r--')
set(gca,'YLim',freqLim)
modPlot
ylabel('cut freq(Hz)')
xlabel('bundle stiff')

figure(f2);
subplot(2,numCols,6)
loglog(cE1.freqs,cE1.sensitivity,'g',cE2.freqs,cE2.sensitivity,'r--');
set(gca,'TickDir','out')
set(gca,'TickLength',[.04 .04])
set(gca,'XLim',[0.01 1000])
set(gca,'XTick',[10.^[-1:1:3]])
set(gca,'YLim',10.^[-6 -3])
set(gca,'YTick',[10.^[-6:1:-3]])

clear cT cE1 cE2 dDai


% Run  'EI kino'
vals    = linspace(0.5e-21,4e-21,50);
exVals  = [1e-21 3.5e-21];
cT      = c;

for j = 1:length(vals)

    cT.EI_kino      = vals(j);

    cT  = numerical_twopart(cT,'torsion spring');
    cT  = calcBunResp(cT,wrt);
    [c_sense,c_cut] = calcPerformance(cT,wrt);

    dDai.val(j)       = vals(j);
    dDai.Dsense(j)    = c_sense;
    dDai.Dcut(j)      = c_cut;

end

cE1             = c;
cE1.EI_kino     = exVals(1);
cE1.freqs       = [10.^linspace(-2,3,500)]';
cE1             = numerical_twopart(cE1,'torsion spring');
cE1             = calcBunResp(cE1,wrt);

cE2             = c;
cE2.EI_kino     = exVals(2);
cE2.freqs       = cE1.freqs;
cE2             = numerical_twopart(cE2,'torsion spring');
cE2             = calcBunResp(cE2,wrt);

figure(f1)
subplot(numRows,numCols,11)
plot(dDai.val,dDai.Dsense); hold on
plot([exVals(1) exVals(1)],[senseLim(1) senseLim(2)],'g--')
plot([exVals(2) exVals(2)],[senseLim(1) senseLim(2)],'r--')
set(gca,'YLim',senseLim)
modPlot
ylabel('max sensitivity (s)')
xlabel('EI kino')
subplot(numRows,numCols,numCols+11);hold on
plot([exVals(1) exVals(1)],[freqLim(1) freqLim(2)],'g--')
plot([exVals(2) exVals(2)],[freqLim(1) freqLim(2)],'r--')
plot(dDai.val,dDai.Dcut)
set(gca,'YLim',freqLim)
modPlot
ylabel('cut freq(Hz)')
xlabel('EI kino')


figure(f2);
subplot(2,numCols,7)
loglog(cE1.freqs,cE1.sensitivity,'g',cE2.freqs,cE2.sensitivity,'r--');
set(gca,'TickDir','out')
set(gca,'TickLength',[.04 .04])
set(gca,'XLim',[0.01 1000])
set(gca,'XTick',[10.^[-1:1:3]])
set(gca,'YLim',10.^[-6 -3])
set(gca,'YTick',[10.^[-6:1:-3]])

clear cT cE1 cE2 dDai


function modPlot
set(gca,'TickDir','out')
set(gca,'TickLength',[.04 .04])
%set(gca,'YLim',[0 max(heights)])
%set(gca,'YTick',[0:1e-5:max(heights)])
tmp=get(gca,'XLim');
set(gca,'XLim',[tmp(1)-range(tmp)/20 tmp(2)]); clear tmp
tmp=get(gca,'YLim');
set(gca,'YLim',[tmp(1)-range(tmp)/20 tmp(2)]); clear tmp
axis square


function c = c_default_theo
%Parameters for all anlayses
c.freqs         = [10.^linspace(-4,3,500)]';
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


function c = c_default_will
%Mean morphometric values from Will's data
c.freqs         = [10.^linspace(-4,3,100)]';
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
c.baseDiameter 	= 9.8e-6;
c.midDiameter 	= 10.9e-6 ;
c.kinoHeight 	= 21.7e-6;
c.cupHeight     = 35.7e-6;
c.numHairs      = 9; %Harris_2003



function [maxSense,cutFreq] = calcPerformance(c,wrt)

maxSense = max(c.sensitivity);

if strcmp(wrt,'freestream')
    
    cutFreq = c.freqs(min(find(c.sensitivity==maxSense)));
    
    if 0
        figure;
        plot(log10(c.freqs),log10(c.sensitivity),'b-',...
            log10(cutFreq),log10(maxSense),'ro')
        grid on
    end
    
elseif strcmp(wrt,'bundle')
    iLow    = find(c.freqs < 0.3);
    iHigh   = find(c.freqs > 50);

    f       = log10(c.freqs);
    s       = log10(c.sensitivity);

    cLow    = polyfit(f(iLow),s(iLow),1);
    cHigh   = polyfit(f(iHigh),s(iHigh),1);

    cutFreq = 10.^((-cLow(2) + cHigh(2)) / (cLow(1) - cHigh(1)));

    if 0
        figure
        plot(f,s,'',f,polyval(cLow,f),'r:',...
            f,polyval(cHigh,f),'r:',...
            log10(cutFreq),polyval(cLow,log10(cutFreq)),'ko')
        grid on
        %pause(.5)
        %close
    end
end

function c = calcBunResp(c,wrt)
% Calculates the frequency response from the data in c for the bundle height.

if ~isfield(c,'M')
    error('need to run numerical solver first')
end

% grab vales from c
h       = c.bunHeight;
M       = c.M;
flwAmp = c.dispAmp;
freqs   = c.freqs;
heights = c.heights;


%Loop thru freqs, calculate performance at bundle
for ii = 1:size(M,2)
    omega       = 2*pi*freqs(ii);
    delta       = (2*c.mu/(omega*c.rho))^0.5;
    
    if strcmp(wrt,'bundle')
        stimFlw   = i*omega*flwAmp * (1-exp(-h*(1+i)/delta));
    elseif strcmp(wrt,'freestream')
        stimFlw   = i*omega*flwAmp;
    end
    
    bundleAmp   = M(1,ii);
    bundleResp  = bundleAmp / stimFlw;
    
    c.sensitivity(ii,1) = abs(bundleResp);
    c.phase(ii,1)       = angle(bundleResp) / pi * 180;
    
    clear omega delta bundleSpd bundleAmp bundleResp
end



function cMap = makeColorMap(numrows,colors)
%makes a colormap given a matrix of that describes a series of colors with rgb values in colums
%numrows specifies the length of the colormap
if ~(numrows/2==floor(numrows/2))
    error('You must have an even number of rows');
end
if size(colors,1)<2
    error('You must give at least 2 colors');
end
reds=[];greens=[];blues=[];
for i = 1:size(colors,1)-1
    reds   = [reds; giveSeries(colors(i,1),colors(i+1,1),round(numrows/(size(colors,1)-1)))];
    greens = [greens; giveSeries(colors(i,2),colors(i+1,2),round(numrows/(size(colors,1)-1)))];
    blues  = [blues; giveSeries(colors(i,3),colors(i+1,3),round(numrows/(size(colors,1)-1)))];
end
cMap = [reds greens blues];
if size(cMap,1)<numrows
    fillOnes    = ones(numrows-size(cMap,1),1);
    cMap        = [cMap; [cMap(end,1)*fillOnes cMap(end,2)*fillOnes cMap(end,3)*fillOnes]];
elseif size(cMap,1)>numrows
    cMap        = cMap(1:numrows,:);
end

function s = giveSeries(num1,num2,numSteps)
if num1==num2
    s   = num1 .* ones(numSteps,1);
else
    s   = [num1 : (num2-num1)/(numSteps-1) : num2]';
end
