function analyze_extract(fname,pathname)



%% Parameter values 

% Command to save a csv file of the data trace
save_csv = 0;

% Command to plot signal data
plot_signal = 0;

% Comparison with Seiste's analysis of Kalmijn trace
kalmijn_compare = 1;

% Viscosity [Pa.s] 
mu =  0.001; 

% Density (kg/m^3)
rho = 1000;	 

% Canal neuromast parameters
CN.ff1 = 10.3;
CN.tau = 0.003;
CN.pc1 = 80;

% Canal radius (m)
CN.R = 0.5e-3; 

% Canal cut-off (Hz)
CN.f_can = 20; 

% Cupular Youngs modulus (Pa)
SN.E = 80; 

% SN cupula radius (m)
SN.r = 4.5e-6; 

% Height of superficial neuromast (m)
SN.h = 50e-6;

% Bundle height (m)
SN.h_bun = 5.2e-6;

% Linear hair bundle stiffness (van Netten & Kroese, 1987) (N/m)
SN.k_hairs     = 0.001;

% Number of hair bundles
SN.num_hairs = 10;


%% Load acquired data

% Prompt for data file
if nargin < 2
    [fname,pathname,findex] = uigetfile('*.mat;','Pick a data file');
end

% Load 'sietse' structure (Siets's data for comparison
if kalmijn_compare
    load(['/Users/mmchenry/Dropbox/Manuscripts/Chapter 5/' ...
          'literature data of sources/sietse_data.mat']);
end

% Load d structure
load([pathname filesep fname])

% Check for dimensions field
if ~isfield(d,'dimen')
    error('Dimension field is undefined, re-run extract_data')
end

clear findex 


%% Translate signal data

% Pad time vector with n values at start
n_vals   = 1;
diff_x   = mean(diff(d.x_vals));
d.x_vals = [min(d.x_vals)+diff_x.*[-n_vals:-1]'; d.x_vals];

% Define time signal vector
s.t = (d.x_vals - min(d.x_vals));

% Define displacement, velocity, acceleration discretely
if strcmp(d.dimen,'Velocity')
    
    % Add zero values at start
    d.y_vals = [zeros(n_vals,1); d.y_vals];
    
    s.vel   = d.y_vals;  
    s.disp  = cumtrapz(s.t,s.vel);
    s.acc   = [0; diff(s.vel)./diff(s.t)]; 
    
elseif strcmp(d.dimen,'Acceleration')
  
    % Add zero values at start
    d.y_vals = [zeros(n_vals,1); d.y_vals];
    
    s.acc   = d.y_vals;
    s.vel   = cumtrapz(s.t,s.acc); 
    s.disp  = cumtrapz(s.t,s.vel);
   
elseif strcmp(d.dimen,'Displacement')
    
    % Add values at start
    d.y_vals = [d.y_vals(1).*ones(n_vals,1); d.y_vals];
    
    s.disp   = d.y_vals;
    s.vel    = [0; diff(s.disp)./diff(s.t)];
    s.acc    = [0; diff(s.vel)./diff(s.t)];  
end

% Trim data to nearest power of 2
two_vals = 2.^[0:20];

% Half the number of samples
s.N = two_vals(find(length(s.t)<two_vals,1,'first')-2);

% Trim vectors to nearest power of 2
s.t    = s.t(1:2*s.N-n_vals);
s.vel  = s.vel(1:2*s.N-n_vals);
s.acc  = s.acc(1:2*s.N-n_vals);
s.disp = s.disp(1:2*s.N-n_vals);

% Add to the end of the time vector
dif_t = mean(diff(s.t));
s.t = [s.t; max(s.t)+dif_t.*[1:n_vals]'];
clear dif_t

% Add values to end of the flow recordings
s.vel  = [s.vel;  zeros(n_vals,1)];
s.acc  = [s.acc;  zeros(n_vals,1)];
s.disp = [s.disp; s.disp(end).*ones(n_vals,1)];

% Spectra of velocity and acceleration
s.V   = fft(s.vel);
s.ACC = fft(s.acc);

% Duration of recording
s.T = range(s.t);

% Inverse of recording duration 
s.df = 1/s.T;

% Index for time series
s.n = [0:(2*s.N-1)]';

% Index for frequency
s.f = [0:(2*s.N-1)]';

% Define units
if strcmp(d.units_y,'au')
    s.unit_disp = 'au';
    s.unit_vel  = 'au';
    s.unit_acc  = 'au';
else
    s.unit_disp = 'm';
    s.unit_vel  = 'm/s';
    s.unit_acc  = 'm/s^2';
end

% Save csv file (for sietse)
if save_csv   
    csvwrite([pathname filesep fname '.csv'],[s.t s.disp s.vel s.acc]);
    disp(' ')
    disp('Saved file:')
    disp([pathname filesep fname '.csv'])
end

clear diff_x findex two_vals t_const


%% Plot signal data

if plot_signal
    figure
    subplot(3,1,1)
    plot(s.t,s.disp,'k')
    %hold on;plot(s.t,cumtrapz(s.t,s.vel),'r--')
    xlabel('Time (s)')
    ylabel(['Displacement' s.unit_disp])
    
    subplot(3,1,2)
    plot(s.t,s.vel,'k')
    %hold on;plot(s.t,[0;diff(s.disp)./diff(s.t)],'r--')
    xlabel('Time (s)')
    ylabel(['Velocity ' s.unit_vel])
    
    subplot(3,1,3)
    plot(s.t,s.acc,'k')
    %hold on;plot(s.t,[0;diff(s.vel)./diff(s.t)],'r--')
    xlabel('Time (s)')
    ylabel(['Acceleration ' s.unit_acc])
end


%% Calc response for CN (wrt flow in canal)

H.df = s.df;
H.f  = s.f;

% Transfer function for CN (wrt flow velocity in canal)
H.CN = (s.f < s.N) .* H_CN(s.f.*s.df,CN.ff1,CN.pc1) + ...
       (s.f > (s.N-1)).* conj(H_CN((2*s.N-s.f).*s.df,CN.ff1,CN.pc1));

% Transfer function for the canal (velocity wrt acceleration)
H.can = (s.f<s.N) .* H_can(s.f.*s.df,CN.f_can) + ...
        (s.f>(s.N-1)) .* conj(H_can((2*s.N-s.f).*s.df,CN.f_can));

% Time signals for response set equal to signal
r.t  = s.t;   
    
% Spectrum of response (hair bundle wrt acceleration)
r.spec_CN = s.ACC .* H.can .* H.CN;

% Response of CN in time domain.
r.resp_CN = real(ifft(r.spec_CN)./s.T);


%% Response for SN


[ra,rphi,ret,H2,rw,rz,rf] = nm4(s.f.*s.df,SN.h_bun,SN.E,SN.r,SN.h,SN.k_hairs,SN.num_hairs);





r.spec_SN = s.V.*H2(:,1);
%spec2 = V.*rf(:,1);

%TODO: calculate appropirate height
%TODO: implement force calculation of response

r.resp_SN = real(ifft(r.spec_SN)./s.T);


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



%% Plot spectra

if 0

% CN Filtering spectra
figure
subplot(2,2,1)
loglog(H.f.*H.df,abs(H.CN),'r--');

subplot(2,2,3)
semilogx(H.f.*H.df,angle(H.CN).*(180/pi),'r--');

subplot(2,2,2)
loglog(H.f.*H.df,abs(H.can),'r--')

subplot(2,2,4)
semilogx(H.f.*H.df,angle(H.can).*(180/pi),'r--');

% Signal spectra
figure
subplot(2,2,1)
loglog(H.f.*H.df,abs(H.resp_CN),'r--');

subplot(2,2,3)
semilogx(H.f.*H.df,angle(H.resp_CN).*(180/pi),'r--');

subplot(2,2,2)
loglog(H.f.*H.df,abs(s.ACC),'r--')

subplot(2,2,4)
semilogx(H.f.*H.df,angle(s.ACC).*(180/pi),'r--');

end


%% Compare Sietse's data with my implementation


if kalmijn_compare
    
    figure
    subplot(3,1,1)
    plot(sietse.t,sietse.vel,'k',s.t,s.vel,'r--')
    title('signal')
    legend('sietse','Mine')
    
    subplot(3,1,2)
    plot(sietse.t,sietse.CN_resp,'k',r.t,r.resp_CN,'r--')
    
    subplot(3,1,3)
    plot(sietse.t,sietse.SN_resp,'k',r.t,r.resp_SN,'r--')
    
    
end

function H_val = H_CN(freq,ft,pc)
% Transfer function describing the freq resp of a CN wrt flow velocity in
% the canal

numer = ( (1/(2*pi.*ft)) + (0.5*sqrt(2)*(1i-1).*sqrt(freq))./(1i*2*pi) ...
        .*(1./ft).^(3/2) - (1/ft).^2 .* freq ./ (3.*1i.*2.*pi));
    
denom = pc + (1i.*freq./ft + 0.5.*sqrt(2).*(1i-1).*(freq./ft).^(3/2) - ...
        ((freq./ft).^2)./3);
    
H_val = numer ./ denom;


function H_out = H_can(freq,f_can)
% Frequency response of the canal.  Defines velocity within the canal wrt
% acceleration at the surface.

H_out = [1./(2.*pi.*f_can.*(1i.*(freq./f_can) + 1))];

 
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








