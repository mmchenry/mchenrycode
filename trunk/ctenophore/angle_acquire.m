function angle_acquire(imPath)


%% Parameters

num_digits = 6;     % Digits at end of image filenames
name_prefix = 'cteno';
t_bin = 1.5;           % Time bin of which to evaluate average data   

%% Define directories

% Prompt for first frame, if not given
if nargin < 1 
    [imPath] = uigetdir(pwd);
end

%% Acquire sequence information

% Look for existing seq_info file and load file
a = dir([imPath filesep 'seq_info.mat']);

if isempty(a)
    %output error message
    error('No sequence data exists for this sequence.');
else
load([imPath filesep 'seq_info']);
end

% Look for existing body_data file and load file
b = dir([imPath filesep 'body_data.mat']);

if isempty(b)
    %output error message
    error('No body data exists for this sequence.');
else
    load([imPath filesep 'body_data']);
end


%% Fill in gaps in the angle data

a = dir([imPath filesep 'plate_data_angl.mat']);

if isempty(a)
    
    % Look for existing plate_data file and load file
    b = dir([imPath filesep 'plate_data.mat']);
    
    if isempty(b)
        %output error message
        error('No plate data exists for this sequence.');
    else
        load([imPath filesep 'plate_data']);
    end
    
    % Step thru each plate
    for i = 1:length(pl)
        
        % Define index where both body and angle data exist
        idx = ~isnan(pl(i).angleY) & ~isnan(body.mouthX) ...
            & ~isnan(body.statX);
        
        % Define index where angle data do not exist, but body data do
        idx_nan = isnan(pl(i).angleY) & ~isnan(body.mouthX) ...
            & ~isnan(body.statX);
        
        % Frame numbers with angle values
        val_frames = find(idx);
        
        % Frame numbers with nans for angle values
        nan_frames = find(idx_nan);
        
        % Cut out nan frames that preceed first value
        nan_frames = nan_frames(find(nan_frames>val_frames(1)));
        
        % Current frame = first frame with angle points
        cFrame = val_frames(1);
        
        % Define current x & y values (global coordinates) for angle coordinate
        xG = pl(i).angleX(cFrame);
        yG = pl(i).angleY(cFrame);
        
        % Define body coordinate system for current frame
        origin = [mean([body.mouthX(cFrame) body.statX(cFrame)]) ...
            mean([body.mouthY(cFrame) body.statY(cFrame)])];
        S = localSystem(origin,[body.mouthX(cFrame) body.mouthY(cFrame)]);
        
        % Define angle coordinate in body system
        [xB,yB] = globalToLocal([xG yG],origin,S);
        
        clear xG yG
        
        % Step through nan frames
        for k = 1:length(nan_frames)
            
            % Define current frame
            cFrame = nan_frames(k);
            
            % Define origin and transform for current frame
            or_c = [mean([body.mouthX(cFrame) body.statX(cFrame)]) ...
                mean([body.mouthY(cFrame) body.statY(cFrame)])];
            S_c  = localSystem(origin,[body.mouthX(cFrame) ...
                body.mouthY(cFrame)]);
            
            % Transform body system point into global coordinates
            [xG,yG] = localToGlobal([xB yB],or_c,S_c);
            
            % Store coordinates
            pl(i).angleX(cFrame) = xG;
            pl(i).angleY(cFrame) = yG;
            
            % If angle is defined in next frame, update coordinate in body
            % coordinates
            if k<length(nan_frames) && idx(cFrame+1)
                
                % Define current x & y values (global coordinates) for angle coordinate
                xT = pl(i).angleX(cFrame+1);
                yT = pl(i).angleY(cFrame+1);
                
                % Define body coordinate system for current frame
                origin = [mean([body.mouthX(cFrame+1) body.statX(cFrame+1)]) ...
                    mean([body.mouthY(cFrame+1) body.statY(cFrame+1)])];
                S = localSystem(origin,[body.mouthX(cFrame+1) body.mouthY(cFrame+1)]);
                
                % Define angle coordinate in body system
                [xB,yB] = globalToLocal([xT yT],origin,S);
                
            end
            
            clear cFrame S_c or_c xT yT xG yG
            
        end
        
        clear xB yB S origin
        
    end
    
    save([imPath filesep 'plate_data_angl'],'pl');
    
else
    clear pl
    disp('Loading plate_data_angl');
    load([imPath filesep 'plate_data_angl'])
    
end


%% Calculate angle data wrt body

a = dir([imPath filesep 'angl_data.mat']);

if isempty(a)
    for i=1:length(pl)
        
        % Define frame numbers with all body and plate data
        fNum = 1:length(pl(i).tipY);
        fNum = fNum(~isnan(pl(i).tipY) & ~isnan(pl(i).baseY) & ~isnan(pl(i).angleY)...
            & ~isnan(body.mouthX) & ~isnan(body.statX));
        
        for j = 1:length(fNum)
            
            % Current frame
            cFrame = fNum(j);
            
            % Find points
            P_tip    = [pl(i).tipX(cFrame)    pl(i).tipY(cFrame)];
            P_base   = [pl(i).baseX(cFrame)   pl(i).baseY(cFrame)];
            P_angl   = [pl(i).angleX(cFrame)  pl(i).angleY(cFrame)];
            
            % Define body coordinate system for current frame
            origin = [mean([body.mouthX(cFrame) body.statX(cFrame)]) ...
                mean([body.mouthY(cFrame) body.statY(cFrame)])];
            S = localSystem(origin,[body.mouthX(cFrame) body.mouthY(cFrame)]);
            
            % Redefine points in body system
            [P_tip(1),P_tip(2)]   = globalToLocal([P_tip(1) P_tip(2)],origin,S);
            [P_base(1),P_base(2)] = globalToLocal([P_base(1) P_base(2)],origin,S);
            [P_angl(1),P_angl(2)] = globalToLocal([P_angl(1) P_angl(2)],origin,S);
            
            % Define postion vectors & angle
            tip_vect  = [-(P_tip(1)-P_base(1))  P_tip(2)-P_base(2)];
            angl_vect = [-(P_angl(1)-P_base(1)) P_angl(2)-P_base(2)];
            
%             angl = (atan2(tip_vect(2),tip_vect(1)) - ...
%                 atan2(angl_vect(2),angl_vect(1))).*180/pi;
            angl = atan2(tip_vect(2),tip_vect(1)) .* 180/pi;
            
            if 0 %(cFrame./str2num(seq.frame_rate)>35) &&  angl>100
                plot([0 tip_vect(1)],[0 tip_vect(2)],'r-',[0 angl_vect(1)],[0 angl_vect(2)],'b');axis equal;grid on
                title(['frame ' num2str(cFrame)  '  ' num2str(angl)])
                ttt=3;
            end
            
            % Store angle
            pl(i).angl(j,1) = angl;
            pl(i).time(j,1) = [cFrame./str2num(seq.frame_rate)];
            
            clear angl T_tip P_base P_angl origin S tip_vect angl_vect
            
        end
        
        clear fNum
    end
    
    save([imPath filesep 'angl_data.mat'],'pl')
    
else
    
    clear pl
    disp('Loading angl_data . . .')
    load([imPath filesep 'angl_data'])
    
end


%% Calculate body angle

% Define index where both body data exist
fNum = 1:length(body.mouthX);
idx = ~isnan(body.mouthX) & ~isnan(body.statX) & ~isnan(pl(1).tipY);

cntr = [mean([body.statX(idx) body.mouthX(idx)],2) ...
    mean([body.statY(idx) body.mouthY(idx)],2)];

% Store timeseries
body.time = [fNum(idx)./str2double(seq.frame_rate)]';
body.angl = (180/pi).*atan2(body.statY(idx)-cntr(:,2),body.statX(idx)-cntr(:,1));
    
body.drive_freq = peakFreq(body.time,body.angl);


clear idx cntr

%% Calculate boxcar averages

% Put together common time intervals

tMin = 0;
tMax = 10^10;

for i = 1:length(pl)   
    tMin = max([tMin pl(i).time(1)]);
    tMax = min([tMax pl(i).time(end)]);   
end

intrvls = tMin:t_bin:tMax;

clear tMin tMax

% Calculate amp and freq over common intervals
for i=1:length(pl)    
    for j = 1:length(intrvls)-1
        idx = (pl(i).time>=intrvls(j)) & (pl(i).time<intrvls(j+1));
        angl = pl(i).angl(idx);
        
        plD(i).t(j)   = mean([intrvls(j) intrvls(j+1)]);
        plD(i).amp(j) = range(angl);
        plD(i).freq(j) = peakFreq(pl(i).time(idx),angl);
    end
    clear time    
end

% Calculate avgerage body angle for same intervals
for j = 1:length(intrvls)-1
    idx = (body.time>=intrvls(j)) & (body.time<intrvls(j+1));
    angl = body.angl(idx);
    
    bD(i).t(j)    = mean([intrvls(j) intrvls(j+1)]);
    bD(i).angl(j) = mean(angl);
end


%% Plot data

figure
xlims = [10^5 0];

subplot(length(pl)+1,1,1)
plot(body.time,body.angl);
grid on
ylabel('Body angl (deg)')

for i=1:length(pl)
    
    subplot(length(pl)+1,1,i+1)
    plot(pl(i).time, pl(i).angl)
    ylabel(['Angl pl' num2str(i)])
    %ylim([0 3])
    grid on
    tmp = xlim;
    %xlims = [min([xlims(1) tmp(1)]) max([xlims(2) tmp(2)])];
    
    clear time 
    %pause
end

figure

xlims = [10^5 0];

subplot(length(pl)+1,1,1)
plot(body.time,body.angl);
grid on
ylabel('Body angl (deg)')
xlabel('time (s)')


for i=1:length(pl)
    
    subplot(length(pl)+1,1,i+1)
    plot(plD(i).t, plD(i).freq,'o-')
    ylabel(['Freq pl' num2str(i)])
    %ylim([0 3])
    grid on
    tmp = xlim;
    %xlims = [min([xlims(1) tmp(1)]) max([xlims(2) tmp(2)])];
    
    clear time 
    %pause
end
xlabel('time (s)')


figure
for i = 1:length(pl)
    h(i) = plot(bD.angl,plD(i).freq,'o');
    set(h(i),'MarkerEdgeColor',give_clr(i))
    set(h(i),'MarkerFaceColor',give_clr(i))
    hold on
    Ltext{i} = ['pl ' num2str(i)];
end
xlabel('Body angle')
ylabel('Beat freq')
xlim([45 135])
set(gca,'XTick',[45 90 135])
axis square
grid on
legend(h,Ltext)

clear h Ltext

figure
for i = 1:length(pl)
    h(i) = plot(bD.angl,plD(i).amp,'o');
    set(h(i),'MarkerEdgeColor',give_clr(i))
    set(h(i),'MarkerFaceColor',give_clr(i))
    hold on
    Ltext{i} = ['pl ' num2str(i)];
end
xlabel('Body angle')
ylabel('Beat amplitude')
axis square
grid on
xlim([45 135])
set(gca,'XTick',[45 90 135])
legend(h,Ltext)

clear h Ltext

return
for i = 1:length(pl)
   subplot(length(pl),1,i)
   %xlim(xlims);
   xlim([0 4.5])
end


function S = localSystem(P1,P2)
% Defines a transformation vector for a local coordinate system in an
% inertial frame of reference.  Uses P1 as the origin and P2 to find the
% direction of the y-axis.  Coordinates must be (1x2) vectors.

if size(P1,1)~=1 || size(P1,2)~=2 ||...
   size(P2,1)~=1 || size(P2,2)~=2
    error('Coordinates must be (1x2) vectors');
end

yAxis       = (P2-P1)./norm(P2-P1);
xAxis       = [yAxis(2); -yAxis(1)];
S           = [xAxis yAxis'];


function [x,y] = localToGlobal(pts,origin,S)

if size(pts,2)~=2 || size(origin,2)~=2 
    error('Coordinates must be a (nx2) vector');
end

pts         = [inv(S)'*pts']';
pts(:,1)    = pts(:,1)+origin(1);
pts(:,2)    = pts(:,2)+origin(2);
x           = pts(:,1);
y           = pts(:,2);


function [x,y] = globalToLocal(pts,origin,S)

if size(pts,2)~=2 || size(origin,2)~=2 
    error('Coordinates must be a (nx2) vector');
end

pts(:,1)    = pts(:,1)-origin(1);
pts(:,2)    = pts(:,2)-origin(2);
pts         = [S'*pts']';
x           = pts(:,1);
y           = pts(:,2);


function [f,P] = fft_data(y,sample_rate)
NFFT = 2^nextpow2(length(y)); % Next power of 2 from length of y
Y = fft(y,NFFT)/length(y);
f = sample_rate/2*linspace(0,1,NFFT/2+1);

P = 2*abs(Y(1:NFFT/2+1));

function fPeak = peakFreq(t,y)

y = (y-mean(y))./range(y);

% Find mean sample rate
dt = mean(diff(t))/10;
samplerate = 1./dt;

t_avg = t(1):dt:t(end);
y_avg = interp1(t,y,t_avg);

% Calculate driving frequency
[f,P] = fft_data(y_avg,samplerate);

fPeak = f(find(P==max(P),1,'first'));


function clr = give_clr(n)
colorset = [...
    0.00 0.00 1.00 ; ... % Data 1 - blue
    0.00 1.00 0.00 ; ...% Data 2 - green
    1.00 0.00 0.00 ; ...% Data 3 - red
    0.00 1.00 1.00 ; ...% Data 4 - cyan
    1.00 0.00 1.00 ; ...% Data 5 - magenta
    0.75 0.75 0.00 ; ...% Data 6 - RGB
    0.25 0.25 0.25 ; ...% Data 7
    0.75 0.25 0.25 ; ...% Data 8
    0.95 0.95 0.00 ; ...% Data 9
    0.25 0.25 0.75 ; ...% Data 10
    0.75 0.75 0.75 ; ...% Data 11
    0.00 0.50 0.00 ; ...% Data 12
    0.76 0.57 0.17 ; ...% Data 13
    0.54 0.63 0.22 ; ...% Data 14
    0.34 0.57 0.92 ; ...% Data 15
    1.00 0.10 0.60 ; ...% Data 16
    0.88 0.75 0.73 ; ...% Data 17
    0.10 0.49 0.47 ; ...% Data 18
    0.66 0.34 0.65 ; ...% Data 19
    0.99 0.41 0.23 ; ...% Data 20
    ];

clr = colorset(n,:);


%figures
% figure
% subplot(5,1,1)
% plot(time(1),angle(1))
% ylabel('Plate 1')
% subplot(5,1,2)
% plot(time(2),angle(2))
% ylabel('Plate 2')
% subplot(5,1,3)
% plot(time(3),angle(3))
% ylabel('Plate 3')
% subplot(5,1,4)
% plot(time(4),angle(4))
% ylabel('Plate 4')
% subplot(5,1,5)
% plot(time(4),angle(5))
% ylabel('Plate 5')
% xlabel('time (s)')
