function play_kalmijn
% Replicates Sietse's mathcad file of Kalmijn's data


%% Parameter values for signal

% Sample period
T = 0.65;

% Number of samples
N = 4096;

% Order of the Bessel function 
bsl_order = 0;

% Decay constant for displacement function
a = 1/700; 

% A constant that factors into the Bessel function
beta = 1;


%% Parameter values for filter

ff1 = 10.3;

tau = 0.003;

pc1 = 80;


%% Calculated parameters & variables

% Samplerate 
df = 1/T;

% Timestep
dt = T/(2*N);


n = [0:(2*N-1)]';

% Frequencies
f = [0:(2*N-1)]';


q = [1:(2*N-1)]';


%% Signal functions

% Components of the displacement waveform
s      = sin(2.8*pi*n./(2*N));
skew   = n.^3.*exp(-a.*n);
g      = besseli(bsl_order,pi*beta.*sqrt(1-(n./N-1).^2))./ ...
        (besseli(bsl_order,pi*beta));

%g = (g-min(g))./range(g);  
    
% Flow displacement
w = s.*skew.*g./(N.^2);

%w = 15.*w./range(w);

% Flow velocity
v_q = [0; [(q>0).*2.*N.*diff(w)./T]];

%v_n = (n>0).*2.*N.*diff(w)./T;


H1 = (f<N).*H(f.*df,ff1,pc1) + (f>(N-1)).* conj(H((2*N-f).*df,ff1,pc1));


% Freq response of the cupula wrt velocity
%H1 = H(f.*df,ff1,pc1);

V = fft(v_q);

% Predicted response in time domain.
spec1 = V.*H1;
resp1 = real(ifft(spec1)./T);




%% Plot results

figure;
subplot(2,1,1)
loglog(f.*df,abs(H1),'b')
%xlim([1 1000])
ylim([.01 10])
set(gca,'YTick',[.01 .1 1 10])
xlabel('Freq (Hz)')
ylabel('Amplitude of response (s)')
title('Freq resp of canal neuromast wrt velocity')

subplot(2,1,2)
h = semilogx(f.*df,angle(H1).*(180/pi));
set(h,'Color','b')
set(gca,'YTick',[-90 -45 0 45])
%xlim([1 1000])
xlabel('Freq (Hz)')
ylabel('Phase of response (deg)')

figure;
subplot(2,1,1)
plot(n.*dt,w,'k')
xlabel('Time (s)')
ylabel('Displacement (m)')

subplot(2,1,2)
plot(n.*dt,v_q./10,'b-',n.*dt,resp1,'r-')
xlabel('Time (s)')
ylabel('Velocity (10 x m/s)')
legend('flow','response')

%ttt= 2;


%% Output excel data

if 0
    out_path = ...
        '/Users/mmchenry/Dropbox/Matlab/van_netten_chapter/numerical_model';
    
    
    header_data = ['n, f, s, skew, g, w, v_q, H1, V, spec1, resp1'];
    out_data    = [n f s skew g w v_q H1 V spec1 resp1];
    
    % Write the column headings (import as csv in excel)
    dlmwrite([out_path filesep 'headers.txt'], ...
        header_data, ...
        'delimiter', '\t')
    
    % Write the data
    csvwrite([out_path filesep 'out data.csv'],out_data);

end


function H_val = H(freq,ft,pc)

numer = 628.*( (1/(2*pi.*ft)) + (0.5*sqrt(2)*(1i-1).*sqrt(freq))./(1i*2*pi) ...
        .*(1./ft).^(3/2) - (1/ft).^2 .* freq ./ (3.*1i.*2.*pi));
    
denom = pc + (1i.*freq./ft + 0.5.*sqrt(2).*(1i-1).*(freq./ft).^(3/2) - ...
        ((freq./ft).^2)./3);
    
H_val = numer ./ denom;