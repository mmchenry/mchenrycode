function numerical_varyParam
% Runs the numerical model of the neuromast


% % Vary EI_kino _________________________________________________________
% figure;
% c       = c_default;
% h       = 5e-6;
% vals    = linspace(0.5e-21,3e-21,4);
% clrs    = makeColorMap(length(vals),[1 0 0;0 0 1]);
% 
% for j = 1:length(vals)
%     c.EI_kino   = vals(j);
%     c           = numerical_twopart(c);   
%     c           = calcFreqResp(c,h);
%     
%     plotFreqResp({c},clrs(j,:))
%     
% end
% 
% subplot(2,1,1)
% title('Varying EI_k_i_n_o');
% clear h vals clrs j c
% 
% % Vary kino height_____________________________________________________
% figure;
% 
% c       = c_default;
% h       = 5e-6;
% vals    = linspace(5e-6,c.cupHeight-1e-6,4);
% clrs    = makeColorMap(length(vals),[1 0 0;0 0 1]);
% 
% for j = 1:length(vals)
%     c.kinoHeight= vals(j);
%     c           = numerical_twopart(c);   
%     c           = calcFreqResp(c,h);
%     
%     plotFreqResp({c},clrs(j,:))
% end
% 
% subplot(2,1,1)
% title('Varying kino height');
% clear h vals clrs j c
% 
% 



% Vary EI tip___________________________________________________________
figure;

c           = c_default;
I_base      = (pi / 64) * c.baseDiameter^4;
c.EI_kino   = 2e-21 +  c.E_matrix .* I_base; 
h           = 5e-6;
vals        = linspace(20,200,4);
clrs        = makeColorMap(length(vals),[1 0 0;0 0 1]);

for j = 1:length(vals)
    c.E_matrix  = vals(j);
    c           = numerical_twopart(c,'torsion spring','special');   
    c           = calcFreqResp(c,h);
    
    plotFreqResp({c},clrs(j,:))
end

subplot(2,1,1)
title('Varying EI_t_i_p');
clear h vals clrs j c






return
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

 disp(' ')
 disp(['This neuromast has a sensitivity of ' num2str(c.peak_sense) ' s'])
 disp(['This peak ocurs at ' num2str(c.peak_freq) ' Hz'])
 disp(' ')

 
 
 
function c = c_default
%Parameters for all anlayses
 c.freqs         = [10.^linspace(-1,3,200)]';
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





function cMap = makeColorMap(numrows,colors)
% Makes a colormap given a matrix of that describes a series of colors with 
% rgb values in colums numrows specifies the length of the colormap

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
