function [heights,amp,phase] = visDeflect(c,freq)
%Visualize the amp and phase of defection for a given frequency.

if ~isfield(c,'M')
    error('need to run numerical solver first')
end

% grab vales from c
bh      = c.bunHeight;
M       = c.M;
freqs   = c.freqs;
heights = c.heights;
L1      = c.kinoHeight;
dispAmp = c.dispAmp;
flwSpd  = i * 2 * pi * freq * dispAmp;

omega       = 2*pi*freq;
delta       = (2*c.mu/(omega*c.rho))^0.5;
stimFlw     = flwSpd * (1-exp(-bh*(1+i)/delta));

% Loop thru heights, evaluate at freq
for ii = 1:size(M,1)
    bDif         = interp1(freqs,M(ii,:),freq);
    bundleResp   = bDif ./ flwSpd;
    amp(ii,1)    = abs(bDif);
    phase(ii,1)  = angle(bundleResp) ./ pi .* 180;

    clear bDif bundleResp
end

%visualize:
if 1
    
    clr = 'b';
    subplot(1,2,1)
    h = plot(10^6 .* amp,10^6 .* heights,'-');hold on
    %set(h,'MarkerSize',2)
    %set(h,'MarkerEdgeColor',clr)
    %set(h,'MarkerFaceColor',clr)
    temp = get(gca,'XTick');
    h2 = plot([min(temp) max(temp)],10^6 .*[L1 L1],'r-');
    ylabel('Height (micrn)');
    xlabel('Amp (micrn)')
    axis square
    hold on
    title(['Stimulus freq. = ' num2str(freq) ' Hz'])
    subplot(1,2,2)
    h = plot(phase,10^6 .* heights,'-');hold on
    %set(h,'MarkerSize',2)
    %set(h,'MarkerEdgeColor',clr)
    %set(h,'MarkerFaceColor',clr)
    temp = get(gca,'XTick');
    h2 = plot([min(temp) max(temp)],10^6 .*[L1 L1],'r-');
    xlabel('Phase (deg)');
    ylabel('Height (micrn)');
    axis square
    hold on
end