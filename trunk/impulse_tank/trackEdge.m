function d = trackEdge(frameRate,calConst)

% Prompt for first tiff in sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[temp,pName,FilterIndex] = uigetfile('*.tif',...
     'Choose first tif image in sequence');
if ~FilterIndex
    error('no file chosen');
end 
[temp,fName,ext,ver] = fileparts([pName temp]);
clear temp

% % Filename for testing code
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fName = 'motor_c001s0002000250';
% ext     = '.tif';
% if ispc
%     pName = '\\atlantis\Users\mjm\Documents\Projects\Impulse_tank\sample_video\motor_C001S0001\';
% else
%     pName = '/Volumes/mjm/Documents/Projects/Impulse_tank/sample_video/motor_C001S0002/';
% end

% Get tiff file sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = getTiffSequence(pName,fName,ext);

% Run calibration, if no constant given
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    [calName,calPath,FilterIndex] = uigetfile('*.tif',...
        'Choose calibration tif image')
    if ~FilterIndex
        error('no file chosen');
    end
    im = imread([calPath calName]);
    calConst = calibrate_image(im)
end
 
% Prompt for frame rate, define time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    temp         = inputdlg('What is this frame rate in fps? ','Frame rate');
    frameRate    = str2num(temp{1});
end
t = linspace(0,length(a)/frameRate,length(a));

% Use datashow to display each frame analyzed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataShow = 0;

% Prompt to get range of x-values to track
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im          = imread([pName a(1).name ext]);
imshow(im);title('Select axis over which to track edge')
[xS,yS]     = choosePoints(im,1,2)
%xS = [104;150];yS = [549;549];
[d1F,xE_start,relPos]   = findPosition(im,xS,yS);
xRange                  = abs(xS(2)-xS(1));
xE                      = xE_start;

if dataShow
    hF          = figure;
    zoomSize    = 100;
    set(hF,'DoubleBuffer','on');
    %set(hF,'Position',[350 90 1515 1025]);
end



for i = 1:length(a)
    im          = imread([pName a(i).name ext]);
    xS          = [round(xE-relPos.*xRange) round(xE+(1-relPos).*xRange)];
    [d1F,xE]    = findPosition(im,xS,yS);
    
    % Display data on video
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if dataShow
        zoomFactor = 15;
        figure(hF)
        subplot(2,1,1)
        imshow(im);hold on
        plot([xS(1) xS(2)],[mean(yS) mean(yS)],'w');
        plot([xE],mean(yS),'+r');
        hold off
        title([num2str(i) ' of ' num2str(length(a)) ' frames'])
        subplot(2,1,2);
        imshow(im);hold on
        plot([xE],mean(yS),'+r');
        hold off
        zoom(zoomFactor)
        xRng    = abs(diff([get(gca,'XLim')]));
        yRng    = abs(diff([get(gca,'YLim')]));
        axis([xE-xRng/2 xE+xRng/2 mean(yS)-yRng/2 mean(yS)+yRng/2]);
        title(['Zoomed ' num2str(zoomFactor) ' x']);
        pause(.2)
    else
        disp([num2str(i) ' of ' num2str(length(a)) ' frames'])
    end
    
    % Store data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d.path          = pName;
    d.calConst      = calConst;
    d.time(i)       = t(i);
    d.frameNum      = a(i).frameNum;
    d.x_pix(i)      = xE;
    d.y_pix(i)      = mean(yS);
    d.x_m(i)        = (xE-xE_start) .* calConst;    
end

figure;
subplot(2,1,1)
plot(1000 .* d.time,1000 .* d.x_m);
grid on;
xlabel('time (ms)');
ylabel('position (mm)');



function a = getTiffSequence(pName,fName_start,ext)
aT          = dir(pName);
fPrefix     = fName_start(1:end-6);
frameNum    = str2num(fName_start(end-5:end));
fIndx       = isfile(fName_start,pName);

if fIndx==0, error('file not found'); end

j = 1;
while 1==1
    %[temp,fName_curr,ext_curr,ver] = fileparts([pName aT.name]);
        temp        = ['000000' num2str(frameNum)];
        tName       = [fPrefix temp(end-5:end)];
        if strcmp(aT(fIndx).name,[tName ext])
            a(j).name       = tName;
            a(j).ext        = ext;
            a(j).path       = pName;
            a(j).frameNum   = frameNum;
            j           = j + 1;
            frameNum    = frameNum + 1;
            fIndx       = fIndx + 1;
            if fIndx > length(aT), break, end
        else 
            break
        end
end




function [d1F,xE,relPos] = findPosition(im,xS,yS)
pixBandWidth    = 2;
edgeSize        = 10; %number of pixels at the edges

%acquire data from images:
x           = [round(min(xS)):round(max(xS))];
y           = [round(mean(yS))-pixBandWidth:round(mean(yS))+pixBandWidth];

d1          = double(  im(y,x)  );

% Average data, add edges:
d1          = mean(d1,1);
d1          = (d1-mean(d1(:))) / range(d1);
x           = [min(x)-edgeSize:max(x)+edgeSize];
d1          = [ones(1,edgeSize).*d1(1) d1 ones(1,edgeSize).*d1(end)];

%Filter & Differentiate
SampleRate  = 1;
cutfreq     = .15;
ff          = cutfreq/(SampleRate/2);
[B A]       = butter(2,ff);
dMean       = mean([d1(:)]);
d1F         = filtfilt(B,A,d1);
d1D			= diff(d1F);

%Normalize derivative:
% base1		= mean(d1D(1:10));
% d1D			= (d1D-base1) ./ range(d1D);
% xD          = x(1:end-1) + mean(diff(x))/2;

%find portion to focus on

%Spline fit
tolerance           = 1.e-1;
%sp          = spaps(x,d1F, tolerance);
sp                  = spmak(augknt(linspace(min(x),max(x),length(d1F)),4),d1F);
dsp                 = fnder(sp);
[minVal,minD]       = fnmin(dsp,[min(x)+edgeSize max(x)-edgeSize]);

spI                 = spmak(augknt(linspace(min(x),max(x),length(d1F)),4),-d1F);
dspI                = fnder(spI);
[maxVal,maxD]       = fnmin(dspI,[min(x)+edgeSize max(x)-edgeSize]);
%dsp = fnder(sp); dspt = fnval(dsp,s); ddspt = fnval(fnder(dsp),s);

% %Fit a 4th order polynomial:
% coef            = polyfit(x,d1F,10);
% coefD           = polyder(coef);
% coefD2          = polyder(coefD);
% Droots      	= roots(coefD);
% D2atRoots       = polyval(coefD2,Droots);
% % deflect         = Droots(find(D2atRoots<0))+ deflectEst;

%remove edges:
x       = x(edgeSize+1:end-edgeSize);
d1      = d1(edgeSize+1:end-edgeSize);
d1F     = d1F(edgeSize+1:end-edgeSize);
%d1D     = d1D(edgeSize+1:end-edgeSize);

xE      = minD;
temp    = abs(x-xE);
relPos  = min(find(temp == min(temp))) / length(x);

if 0
    figure;
    xTemp = linspace(min(x),max(x),1000);
    subplot(2,1,1);
    plot(x,d1,'.k');grid on;hold on
    plot(x,d1F,'b-')
    plot(xTemp,fnval(sp,xTemp),'g-')
    plot([xE xE],[min(d1) max(d1)],'r-')
    subplot(2,1,2);
    %plot(xD,d1D,'b.');grid on;hold on
    plot(xTemp,fnval(dsp,xTemp),'r-');grid on;hold on
    plot([xE xE],[-1 1],'r-')
    xlabel('x position (pix)');
end



function  indx = isfile(fName,fPath)
a	= dir(fPath);
y	= 0;
for i = 1:length(a)
    if (length(a(i).name) > 3) & (a(i).name(end-3)=='.')
        cName	= a(i).name(1:end-4);
    else
        cName   = a(i).name;
    end
    if (length(a(i).name) > 3) & strcmp(cName,fName)
        y = 1;
        break
    end
end
if ~y
    indx = 0;
else
    indx = i;
end


function [x,y] = choosePoints(img,link,numPts)
%Used for finding coordinate points on a static image 'img'.
figure;
imshow(img);
hold on;
set(gcf,'DoubleBuffer','on');
disp(' '); disp(' ');
disp('Left mouse button picks points.');disp(' ');
disp('Right mouse button removes last point.');disp(' ');
disp('Press return to stop.')
n = 0;
but = 1;
while 1 == 1
    [xi,yi,but] = ginput(1);
    imshow(img);hold on
    if isempty(but)
        break
    elseif but==1 %Left click
        n = n+1;
        if n > numPts
            n = numPts;
        end
        x(n) = xi;
        y(n) = yi;
        if link
            plot(x,y,'ro-')
        else
            plot(x,y,'ro')
        end
    elseif but==3 %Right click
        if n-1 < 1
            n = 0;
            x = [];
            y = [];
        else
            n = n-1;
            x = x(1:n);
            y = y(1:n);
        end
        if link
            plot(x,y,'ro-')
        else
            plot(x,y,'ro')
        end
    end
    hold off
end
x = x'; y = y';
if ~(length(x)==numPts)
    disp(' ');disp(' '); warning(['Try again -- you need to find ' num2str(numPts) ' points']);
    [x,y] = choosePoints(img,link,numPts);
end
close