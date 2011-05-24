function spir_analysis(dPath)


%% Set paths

if nargin < 1
    %dPath = 'C:\Users\Guest\Documents\spir\recordings\C001-R001-17-Jul-2010';
    dPath = 'C:\Users\Guest\Documents\spir\recordings\R005-sp003-r002-29-Jul-2010';
end

%% Find out desired mode

but = questdlg('What kind of experiment do you want to analyze?',...
               'Experiment type','2 camera','camera + sensor','2 camera');

switch but
%% 2 Camera experiment
case '2 camera'
    
    % Set frame rate
    frame_rate = 15;
    
    % Check for directories
    a0 = dir([dPath filesep 'pike']);
    a1 = dir([dPath filesep 'marlin']);
    if isempty(a0)
        error(['There needs to be a directory called "pike" within ' ...
               'the root directory'])
    elseif isempty(a1)
        error(['There needs to be a directory called "marlin" within ' ...
               'the root directory'])
    end
    
    % Load pixel data from cam0
    load([dPath filesep 'pike' filesep 'pixel_data.mat'])
    d0 = d;
    clear d
    
    % Load pixel data from cam1
    load([dPath filesep 'marlin' filesep 'pixel_data.mat'])
    d1 = d;
    clear d
    
    % Create time vectors
    t0 = 0:1/frame_rate:(length(d0.pixVal)-1)/frame_rate;
    t1 = 0:1/frame_rate:(length(d1.pixVal)-1)/frame_rate;
    
    figure;
    plot(t0,(d0.pixVal-mean(d0.pixVal))./range(d0.pixVal),'+-')
    hold on
    plot(t0,(d1.pixVal-mean(d1.pixVal))./range(d1.pixVal),'ro-')
    grid on
    hold off
    xlabel('Time (s)')
    ylabel('Normalized pixel intensity')
    legend('cam0','cam1')

    
%% Camera + sensor experiment
case 'camera + sensor'               

% Set parameters

% Sample rate for daq (Hz)
%sample_rate = 1000;

% Frame rate for video (fps)
%frame_rate  = 15;

cut_freq_lowP = 1;

cut_freq_highP = 10;


% Load data

% load 'data' vector (from daq)
warning off
load([dPath filesep 'recording_data.mat'])
warning on

data = s.ch_data;
sample_rate = s.sample_rate;
frame_rate = s.frame_rate;

clear ai0 chan s

t_daq = (0:1/sample_rate:(length(data)-1)/sample_rate);

data_f = butter_filt(data,sample_rate,cut_freq_lowP,'low');  

data_high = butter_filt(data,sample_rate,cut_freq_highP,'high'); 

% Load video
load([dPath filesep 'pixel_data.mat'])

t_vid = 0:1/frame_rate:(length(d.pixVal)-1)/frame_rate;

% Load video data
%load([dPath filesep 

figure;
% subplot(4,1,1)
% plot(t_daq,data)
% hold on
% xlabel('Time (s)')
% ylabel('V')
% 
% 
% subplot(4,1,1)
% plot(t_daq,data_f,'r-')
% xlabel('Time (s)')
% ylabel('V')
% 
% subplot(4,1,2)
% fft_plot(data_high,sample_rate);
% 
% subplot(4,1,3)
% plot(t_vid,d.pixVal)
% 
% subplot(4,1,4)
subplot(2,1,1)
plot(t_daq,(data_f-mean(data_f))./range(data_f),'r-')
hold on
plot(t_vid,(d.pixVal-mean(d.pixVal))./range(d.pixVal),'-')
xlabel('Time (s)')
ylabel('Normalized intensity')
hold off

[fv,Pv] = fft_data((d.pixVal-mean(d.pixVal))./range(d.pixVal),frame_rate);
[fd,Pd] = fft_data((data_f-mean(data_f))./range(data_f),sample_rate);

% Plot single-sided amplitude spectrum.
subplot(2,1,2)
plot(fv,Pv,'b')
hold on
plot(fd,Pd,'r-')
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
xlim([0 3])

legend('Video','Sensor')

end % switch

function data_filtered = butterworth(data,sample_rate,cutfreq,type)
%Returns data filtered by a butterworth filter
%  data - vector of data for filtering
%  sample_rate - rate of sampling for data (must have equivalent intervals)
%  cutfreq - cut off frequency (must be 2 values for bandstop)
%  type - 'low' for lowpass (default), 'high' for highpass, 'stop' for bandstop

if nargin < 4
    type = 'low';
end

if strcmp(type,'stop') && ~(length(cutfreq)==2)
    error('Stop pass must have a two-element cutfreq')
end 

ff              = cutfreq./(sample_rate./2);
[B A]           = butter(2,ff,type);
data_filtered   = filtfilt(B,A,data);     % Filtered data


function data_filtered = butter_filt(data,sample_rate,cut_freq,type) 
% High-pass or low-pass butterworth filter

% All frequency values are in Hz.

% Nyquist freq.
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


function [f,P] = fft_data(y,sample_rate)
NFFT = 2^nextpow2(length(y)); % Next power of 2 from length of y
Y = fft(y,NFFT)/length(y);
f = sample_rate/2*linspace(0,1,NFFT/2+1);

P = 2*abs(Y(1:NFFT/2+1));

