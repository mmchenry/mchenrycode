function play_fft



Fs = 125;    %sample rate
Fo = 0.5;      % Freq of oscillations
lag = pi*(90)/180;   % Lag of second curve wrt first

time = 0:(1/Fs):7;
L = length(time);

y1 = 5.*cos(2.*pi.*time*Fo);
y2 = 5.*cos((2.*pi).*time*Fo-lag);



% APPROACH 1: Does better on phase

Y1 = fft(y1)/L;
Y2 = fft(y2)/L;
len = floor(length(Y1)/2)
f = Fs/2*linspace(0,1,len);

Y1_amp = 2*abs(Y1(1:len));
Y2_amp = 2*abs(Y2(1:len));

% Y1_ang = 180.*unwrap(angle(Y1(1:NFFT/2)))./pi;
% Y2_ang = 180.*unwrap(angle(Y2(1:NFFT/2)))./pi;

Y1_ang = 180.*(angle(Y1(1:len)))./pi;
Y2_ang = 180.*(angle(Y2(1:len)))./pi;

phs = 180.*(angle(Y1(1:len)./Y2(1:len)))./pi;

figure;
subplot(3,1,1)
plot(time,y1,time,y2,'r-')
grid on

subplot(3,1,2)
plot(f,Y1_amp,'b',f,Y2_amp,'r-')
set(gca,'XLim',[0 20])
grid on

subplot(3,1,3)
%plot(f,Y1_ang,'b-+',f,Y2_ang,'r-+')
plot(f,phs,'k+-')
set(gca,'XLim',[0 20])
grid on


% APPROACH 2: Does better on amplitude

NFFT = 2^nextpow2(L)
Y1 = fft(y1,NFFT)/L;
Y2 = fft(y2,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2);

Y1_amp = 2*abs(Y1(1:NFFT/2));
Y2_amp = 2*abs(Y2(1:NFFT/2));

Y1_ang = 180.*(angle(Y1(1:NFFT/2)))./pi;
Y2_ang = 180.*(angle(Y2(1:NFFT/2)))./pi;

figure;
subplot(3,1,1)
plot(time,y1,time,y2,'r-')
grid on

subplot(3,1,2)
plot(f,Y1_amp,'b',f,Y2_amp,'r-')
set(gca,'XLim',[0 20])
grid on

subplot(3,1,3)
plot(f,Y1_ang,'b',f,Y2_ang,'r-')
set(gca,'XLim',[0 20])
grid on


return
% Example from matlab demos

t = 0:1/100:10-1/100;
x = sin(2*pi*15*t) + sin(2*pi*40*t);
L = length(t);

NFFT = 2^nextpow2(L)
y = fft(x,NFFT); 
m = abs(y);
p = unwrap(angle(y));

f = (0:length(y)-1)'*100/length(y);
subplot(3,1,1), plot(t,x)
subplot(3,1,2), plot(f,m), 
ylabel('Abs. Magnitude'), grid on
subplot(3,1,3), plot(f,p*180/pi)
ylabel('Phase [Degrees]'), grid on
xlabel('Frequency [Hertz]')