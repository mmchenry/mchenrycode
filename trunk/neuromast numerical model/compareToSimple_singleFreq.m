function compareToSimple_singleFreq(c)

if ~isfield(c,'M')
    error('need to run numerical solver first')
end
if size(c.M,2)>1
    error('just one freq');
end

% grab vales from c
bh      = c.bunHeight;
M       = c.M;
dispAmp = c.dispAmp;
freqs   = c.freqs;
heights = c.heights;
mu      = c.mu;
rho     = c.rho;
L       = c.cupHeight;

bh=L;


a       = c.baseDiameter/2;
EI      = c.E_matrix .* (pi*a^4/4);
omega   = 2 * pi * freqs;
bDif    = interp1(heights,M,bh);
flwSpd  = i * 2 * pi * freqs * dispAmp;
sense   = abs(bDif) / abs(flwSpd)
phase   = (angle(bDif)) / pi * 180
cMoment = real( c.ddM .* EI );


%simple:
x       = L/2;
kn      = c.bundleStiff .* c.numHairs;
D0      = dispAmp;
delta   = ( 2 * mu / rho / omega )^0.5;
s       = a ./ 2^.5 / delta;
g       = 0.577 + log(s);
G       = -g / (g^2 + (pi/4)^2);
d       = D0 * (1 + i) * bh ./ delta;


%version for negative intensity equation:
%q       = -4 * pi * mu * G * d * i *omega;
%moment  = real( ( (2/3 - 2*i/3) *D0*G*pi .* (L-x).^2 .* (2*L+x) .* mu*omega) ./ ...
%            delta);
% v  = ( (1/30 - i/30) .* D0.*G.*pi.*x .* (40*EI.*x.^3 + kn.*x.* ...
%      (20*L^3 - 10*L^2.*x + x.^3)) .*mu*omega) ./ (EI*kn*delta); %version with neg. intensity 

moment  = - ( (2/3 - 2*i/3) *D0*G*pi * (L-x)^2 * (2*L+x) * mu*omega) / ...
            delta;
        
v  = -( (1/30 - i/30) * D0*G*pi.*x.^2 * (20*L^3 - 10*L^2.*x + x.^3) .* mu*omega )./ ...
    ( EI*delta);

sense_v = abs(v) / abs(flwSpd)
phase_v = (angle(v)) / pi * 180

%  figure;
%  plot(heights,cMoment,'b',heights,moment,'r-')