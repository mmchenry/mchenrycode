function d = sensitivityAnalysis(wrt)
%Uses TWO-PART, TORSION SPRING, NUMERICAL

if nargin < 1
    wrt = 'bundle';
end

if ~( strcmp(wrt,'bundle') | strcmp(wrt,'freestream'))
    error('wrt input not an option');
end

c   = c_default_theo;
jj  = 1;
fracChange = .05;

if strcmp(wrt,'freestream')
    c.freqs         = [10.^linspace(-0.5,1.5,500)]';
end

% Vary baseDiameter _________________________________________________

cT  = c;
cT.baseDiameter = c.baseDiameter.*(1-fracChange);

cL  = numerical_twopart(cT,'torsion spring');
cL  = calcBunResp(cL,wrt);
[cL_sense,cL_cut] = calcPerformance(cL,wrt);


cT  = c;
cT.baseDiameter = c.baseDiameter .* (1+fracChange);

cH  = numerical_twopart(cT,'torsion spring');
cH  = calcBunResp(cH,wrt);
[cH_sense,cH_cut] = calcPerformance(cH,wrt);


d.params{jj} = 'base diameter';
d.Dsense(jj) = (cH_sense - cL_sense) / cL_sense;
d.Dcut(jj) = (cH_cut - cL_cut) / cL_cut;

jj = jj+1;
clear cL cH


% Vary midDiameter __________________________________________________

cT  = c;
cT.midDiameter = c.midDiameter.*(1-fracChange);

cL  = numerical_twopart(cT,'torsion spring');
cL  = calcBunResp(cL,wrt);
[cL_sense,cL_cut] = calcPerformance(cL,wrt);


cT  = c;
cT.midDiameter = c.midDiameter .* (1+fracChange);

cH  = numerical_twopart(cT,'torsion spring');
cH  = calcBunResp(cH,wrt);
[cH_sense,cH_cut] = calcPerformance(cH,wrt);


d.params{jj} = 'mid diameter';
d.Dsense(jj) = (cH_sense - cL_sense) / cL_sense;
d.Dcut(jj) = (cH_cut - cL_cut) / cL_cut;

jj = jj+1;
clear cL cH


% Vary numHairs ____________________________________________________

cT  = c;
cT.numHairs = c.numHairs - 1;

cL  = numerical_twopart(cT,'torsion spring');
cL  = calcBunResp(cL,wrt);
[cL_sense,cL_cut] = calcPerformance(cL,wrt);


cT  = c;
cT.numHairs = c.numHairs + 1;

cH  = numerical_twopart(cT,'torsion spring');
cH  = calcBunResp(cH,wrt);
[cH_sense,cH_cut] = calcPerformance(cH,wrt);


d.params{jj} = 'num hairs';
d.Dsense(jj) = (cH_sense - cL_sense) / cL_sense;
d.Dcut(jj) = (cH_cut - cL_cut) / cL_cut;

jj = jj+1;
clear cL cH


% Vary kinoHeight _________________________________________________

cT  = c;
cT.kinoHeight = c.kinoHeight .* (1-fracChange);

cL  = numerical_twopart(cT,'torsion spring');
cL  = calcBunResp(cL,wrt);
[cL_sense,cL_cut] = calcPerformance(cL,wrt);


cT  = c;
cT.kinoHeight = c.kinoHeight .* (1+fracChange);

cH  = numerical_twopart(cT,'torsion spring');
cH  = calcBunResp(cH,wrt);
[cH_sense,cH_cut] = calcPerformance(cH,wrt);


d.params{jj} = 'kino height';
d.Dsense(jj) = (cH_sense - cL_sense) / cL_sense;
d.Dcut(jj) = (cH_cut - cL_cut) / cL_cut;

jj = jj+1;
clear cL cH


% Vary cupHeight __________________________________________________

cT  = c;
cT.cupHeight = c.cupHeight .* (1-fracChange);

cL  = numerical_twopart(cT,'torsion spring');
cL  = calcBunResp(cL,wrt);
[cL_sense,cL_cut] = calcPerformance(cL,wrt);


cT  = c;
cT.cupHeight = c.cupHeight .* (1+fracChange);

cH  = numerical_twopart(cT,'torsion spring');
cH  = calcBunResp(cH,wrt);
[cH_sense,cH_cut] = calcPerformance(cH,wrt);


d.params{jj} = 'cupula height';
d.Dsense(jj) = (cH_sense - cL_sense) / cL_sense;
d.Dcut(jj) = (cH_cut - cL_cut) / cL_cut;

jj = jj+1;
clear cL cH


% Vary E_matrix __________________________________________________

cT  = c;
cT.E_matrix = c.E_matrix .* (1-fracChange);

cL  = numerical_twopart(cT,'torsion spring');
cL  = calcBunResp(cL,wrt);
[cL_sense,cL_cut] = calcPerformance(cL,wrt);


cT  = c;
cT.E_matrix = c.E_matrix .* (1+fracChange);

cH  = numerical_twopart(cT,'torsion spring');
cH  = calcBunResp(cH,wrt);
[cH_sense,cH_cut] = calcPerformance(cH,wrt);


d.params{jj} = 'E matrix';
d.Dsense(jj) = (cH_sense - cL_sense) / cL_sense;
d.Dcut(jj) = (cH_cut - cL_cut) / cL_cut;

jj = jj+1;
clear cL cH


% Vary EI_kino __________________________________________________

cT  = c;
cT.EI_kino = c.EI_kino .* (1-fracChange);

cL  = numerical_twopart(cT,'torsion spring');
cL  = calcBunResp(cL,wrt);
[cL_sense,cL_cut] = calcPerformance(cL,wrt);


cT  = c;
cT.EI_kino = c.EI_kino .* (1+fracChange);

cH  = numerical_twopart(cT,'torsion spring');
cH  = calcBunResp(cH,wrt);
[cH_sense,cH_cut] = calcPerformance(cH,wrt);


d.params{jj} = 'EI kino';
d.Dsense(jj) = (cH_sense - cL_sense) / cL_sense;
d.Dcut(jj) = (cH_cut - cL_cut) / cL_cut;

jj = jj+1;
clear cL cH


% Visualizing results [bars]_________________________________

figure;
subplot(2,1,1)
bar(d.Dsense)
ylabel('sensitivity')
subplot(2,1,2)
bar(d.Dcut)
ylabel('cut freq')




% Visualizing results [color fields]_______________________________
if 0
clrs = [.3 0 .4; 1 1 1; 0 .2 0];
cMap = makeColorMap(128,clrs);
rng  = 0.7;

%figure
for ii = 1:length(d.params)
    clr = cMap(round(size(cMap,1) * (d.Dsense(ii)+rng)/(2*rng)),:);
    d.clr_Dsense(ii,:) = clr;
%     subplot(length(d.params),1,ii)
%     h = fill([0 0 1 1],[0 1 1 0],clr);
%     set(h,'EdgeColor','none')
%     axis square
%     set(gca,'Visible','off')
    %ylabel(d.params{ii});
end

%figure
for ii = 1:length(d.params)
    clr = cMap(round(size(cMap,1) * (d.Dcut(ii)+rng)/(2*rng)),:);
    d.clr_Dcut(ii,:) = clr;
%     subplot(length(d.params),1,ii)
%     h = fill([0 0 1 1],[0 1 1 0],clr);
%     set(h,'EdgeColor','none')
%     axis square
%     set(gca,'Visible','off')
    %ylabel(d.params{ii});
end

figure;
h = colorbar('horiz'); 
set(h,'EdgeColor','none')
set(h,'TickLength',[0 0])
set(gca,'Visible','off')
colormap(makeColorMap(64,clrs))

end



function c = c_default_theo
%Parameters for all anlayses
c.freqs         = [10.^linspace(-2,3,500)]';
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
