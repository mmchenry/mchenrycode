function compareToSimple(c)



if ~isfield(c,'M')
    error('need to run numerical solver first')
end

% grab vales from c
bh      = c.bunHeight;
%bh      = c.kinoHeight;
M       = c.M;
dispAmp = c.dispAmp;
freqs   = c.freqs;
heights = c.heights;

EI      = c.E_matrix .* (pi*(c.baseDiameter/2)^4/4);

%Loop thru freqs, calculate performance at bundle
for ii = 1:size(M,2)
    omega       = 2 * pi * freqs(ii);
    v(ii,1)     = simpleBeam(bh,omega,c.mu,c.rho,c.cupHeight,...
                    dispAmp,c.baseDiameter/2,EI);
    bDif(ii,1)   = interp1(heights,M(:,ii),bh);
    flwSpd      = i * 2 * pi * freqs(ii) * dispAmp;
    %flwSpd      =  freqs(ii).*dispAmp;
    sense(ii,1)  = abs(bDif(ii)) / abs(flwSpd);
    sense_v(ii,1)= abs(v(ii)) / abs(flwSpd);
    %phase(ii,1)  = (angle(bDif(ii))) / pi * 180;
    %phase(ii,1)  = (angle(bDif(ii)) - angle(flwSpd)) / pi * 180;
    phase(ii,1)  = (angle(bDif(ii))) / pi * 180;
    phase_v(ii,1)= (angle(v(ii))) / pi * 180;
end

%calc sensitivity
c.sensitivity   = max(sense);
c.peak_freq     = freqs(min(find(sense==c.sensitivity)));

%calc cutoff
%Dphase      = 

%visualize:
if 1
    figure;
    clr = 'b';
    subplot(2,2,1)
    h = plot(log10(freqs),log10(sense),'o-'); hold on
    h2 = plot(log10(c.peak_freq),log10(c.sensitivity),'or');
    set(h,'MarkerSize',3)
    set(h,'MarkerEdgeColor',clr)
    set(h,'MarkerFaceColor',clr)
    ylabel('Log10 (Sensivity (s))')
    xlabel('Log10 Freq. (Hz)');
    axis square
    grid on
    subplot(2,2,3)
    h = plot(log10(freqs),phase,'o-');
    set(h,'MarkerSize',3)
    set(h,'MarkerEdgeColor',clr)
    set(h,'MarkerFaceColor',clr)
    ylabel('Phase (deg)');
    xlabel('Log10 Freq. (Hz)');
    axis square
    grid on
    
    
    clr = 'b';
    subplot(2,2,2)
    h = plot(log10(freqs),log10(sense_v),'o-'); hold on
    set(h,'MarkerSize',3)
    set(h,'MarkerEdgeColor',clr)
    set(h,'MarkerFaceColor',clr)
    ylabel('Log10 (Sensivity (s))')
    xlabel('Log10 Freq. (Hz)');
    axis square
    grid on
    subplot(2,2,4)
    h = plot(log10(freqs),phase_v,'o-');
    set(h,'MarkerSize',3)
    set(h,'MarkerEdgeColor',clr)
    set(h,'MarkerFaceColor',clr)
    ylabel('Phase (deg)');
    xlabel('Log10 Freq. (Hz)');
    axis square
    grid on
end




function v = simpleBeam(x,omega,mu,rho,L,D0,a,EI)


delta   = (2 * mu / rho / omega);
s       = a ./ 2^.5 / delta;
g       = 0.577 + log(s);
G       = -g / (g^2 + (pi/4)^2);
d       = D0 * (1 + i) * x ./ delta;
%q       = 4 * pi * mu * G * d * i *omega;
% M       = - ( (2/3 - 2*i/3) *D0*G*pi * (L-y)^2 * (2*L+y) * mu*omega / ...
%             delta;

% v  = -( (1/30 - i/30) * D0*G*pi*x^2 * (20*L^3 - 10*L^2*x + x^3) * mu*omega )/ ...
%     ( EI*delta);
        
v  = -( (1/30 - i/30) * D0*G*pi.*x.^2 * (20*L^3 - 10*L^2.*x + x.^3) .* mu*omega )./ ...
    ( EI*delta);