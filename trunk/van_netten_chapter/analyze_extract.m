function analyze_extract



%% Load acquired data

% Prompt for data file
[fname,pathname,findex] = uigetfile('*.mat;','Pick a data file');

% Load 'sietse' structure
load('/Users/mmchenry/Dropbox/Manuscripts/Chapter 5/literature data of sources/sietse_data.mat')


% Load d structure
load([pathname filesep fname])



%% Parameter values for canal neuromast filter

ff1 = 10.3;
tau = 0.003;
pc1 = 80;


%these are the constant parameters
mu =  0.001; %0.001002;  %viscosity [Pa.s]
rho = 1000;	 %density [kg/m3]
R = 0.5e-3; % Canal radius (m)
f_can = 20; %Hz


%% Parameter values for superficial neuromast filter

% Youngs modulus
E = 80; % Pa

% Cupula radius
r = 4.5e-6; % m

% Height of superficial neuromast
h = 50e-6; % m


k_hairs     = 0.001; %N/m (van Netten & Kroese, 1987) 

num_hairs = 10;

%% Calculated parameters & variables

% Adjust for units for time 
if strcmp(d.units_x,'s')
    t_const = 1;
    
elseif strcmp(d.units_x,'ms')
    t_const = 1./1000;
    
else
    error(['Units of ' d.units_x ' unrecognized'])

end


% Trim data to nearest power of 2
two_vals = 2.^[0:20];

% Half the number of samples
N = two_vals(find(length(d.x_vals)<two_vals,1,'first')-2);


% Sample durtaion
%T = mean(diff(d.x_vals(1:(2*N)))) * t_const;
T = range(d.x_vals(1:(2*N)));

% Samplerate 
df = 1/T;

% Index for time series
n = [0:(2*N-1)]';

% Index for frequency
f = [0:(2*N-1)]';

% Time vector
t = n.*T;

%% Calculate signal properties 


%v_q = d.y_vals(1:(2*N));

v_q = sietse.accel;

if size(v_q,1)<size(v_q,2)
    v_q = v_q';
end

% Freq response of the cupula wrt velocity
%H1 = H(f.*df,ff1,pc1);

% Specturm of signal
V = fft(v_q);

%TODO:ACC


%% Calc response for CN (wrt flow in canal)

% Transfer function for CN (wrt flow velocity in canal)
H1 = (f<N).*H(f.*df,ff1,pc1) + (f>(N-1)).* conj(H((2*N-f).*df,ff1,pc1));

% Transfer function for the canal
delta = sqrt(mu./(rho.*pi.*f.*df));

r = 0;
tmp1 = sqrt(-2*1i) .* r./delta;
tmp2 = sqrt(-2*1i) .* R./delta;
%H_can = (-1i)./(2*pi*f.*df) .* (1 - besselj(0,tmp1)./besselj(0,tmp2));

H3 = (f<N).*H_can(f.*df,f_can) + (f>(N-1)).* conj(H_can((2*N-f).*df,f_can));

% Spectrum of response
spec1 = V.*H3.*H1;

% Response of CN in time domain.
resp1 = real(ifft(spec1)./T);



%% Response for SN


[ra,rphi,ret,H2,rw,rz,rf] = nm4(f.*df,5.2e-6,E,r,h,k_hairs,num_hairs);





spec2 = V.*H2(:,1);
%spec2 = V.*rf(:,1);

%TODO: calculate appropirate height
%TODO: implement force calculation of response

resp2 = real(ifft(spec2)./T);

%% Plot results

if 0
figure;
subplot(4,1,1:2)
plot(t,v_q,'k-')
xlabel('Time (s)')
ylabel(['(' d.units_y ')'])
title(['Time signal for ' fname])

subplot(4,1,3)
loglog(f.*df,abs(V),'k')
xlabel('Freq (Hz)')
ylabel('Amplitude')

subplot(4,1,4)
semilogx(f.*df,(180/pi).*angle(V),'k')
xlabel('Freq (Hz)')
ylabel('Phase (deg)')


figure;
subplot(2,1,1)
h = plotyy(t,v_q,t,resp1);
set(h(1),'YColor','k')
set(h(2),'YColor','r')
set(get(h(1),'YLabel'),'String','Measurement')
set(get(h(2),'YLabel'),'String','Response')
set(get(h(1),'Children'),'Color','k')
set(get(h(2),'Children'),'Color','r')
xlabel('Time (s)')
ylabel(['(' d.units_y ')'])
%title('Signal')
grid on
%legend('measurement','response')


subplot(2,1,2)
h = plotyy(t,v_q,t,resp2);
set(h(1),'YColor','k')
set(h(2),'YColor','r')
set(get(h(1),'YLabel'),'String','Measurement')
set(get(h(2),'YLabel'),'String','Response')
set(get(h(1),'Children'),'Color','k')
set(get(h(2),'Children'),'Color','r')
xlabel('Time (s)')
ylabel(['(' d.units_y ')'])
%title('Signal')
grid on
%legend('measurement','response')

% subplot(2,1,2)
% plot(t,resp1,'r-')
% grid on
% title('Response')
% xlabel('Time (s)')


end




% Filter the data
%y_pos_f = butter_filt(x_vals,y_pos,cut_freq,'low');


%% Compare Sietse's data with my implementation


figure
subplot(3,1,1)
plot(sietse.t,sietse.vel,'k',t,v_q,'r--')
title('signal')
legend('sietse','me')

subplot(3,1,2)
plot(sietse.t,sietse.CN_resp,'k')

%subplot(3,1,3)
hold on
plot(t,resp1,'r--')


ttt=4

 
function data_filtered = butter_filt(time,data,cut_freq,type) 
% High-pass or low-pass butterworth filter
% All frequency values are in Hz.
 
% Check for constant sample rate
if min(round(diff(time)*1000))~=max(round(diff(time)*1000))
    error('Cannot analyze variable sample rate data')
end
 
% Calc sample rate & Nyquist freq.
sample_rate = 1./mean(diff(time));
Nqst = sample_rate/2;   
 
% Calculate stopband frequency
if strcmp(type,'high')
    stop_freq = max([(cut_freq - Nqst/10) .01]);  
 
elseif strcmp(type,'low')
    stop_freq = min([(cut_freq + Nqst/10) (Nqst-.01)]); 
 
end
 
% Stopband Attenuation (dB)
Astop = 10;
 
% Passband Ripple (dB)
Apass = 1;   
 
% Normalise the cutoff freq. to the Nyquist freq
Wp = cut_freq/Nqst;
 
% Normalise the stoppass freq. to the Nyquist freq
Ws = stop_freq/Nqst;
 
% Check cutoff
if (Wp > 1) || (Ws > 1)
    error('Cutoff or bandpass frequencies cannot exceed the Nyquist freq')
end
 
% Calculate the order from the parameters using BUTTORD.
[N,Fc] = buttord(Wp, Ws, Apass, Astop);    
    
% Calculate the zpk values using the BUTTER function.
[B A] = butter(N, Fc, type);
 
% Plot frequency reponse
%freqz(B,A,512,sample_rate); 
 
% Filter the data
data_filtered   = filtfilt(B,A,data); 


function H_val = H(freq,ft,pc)
% Transfer function describing the freq resp of a CN wrt flow velocity in
% the canal

numer = ( (1/(2*pi.*ft)) + (0.5*sqrt(2)*(1i-1).*sqrt(freq))./(1i*2*pi) ...
        .*(1./ft).^(3/2) - (1/ft).^2 .* freq ./ (3.*1i.*2.*pi));
    
denom = pc + (1i.*freq./ft + 0.5.*sqrt(2).*(1i-1).*(freq./ft).^(3/2) - ...
        ((freq./ft).^2)./3);
    
H_val = numer ./ denom;


function H_out = H_can(freq,f_can)

H_out = [1./(2.*pi.*f_can.*(1i.*(freq./f_can) + 1))];





