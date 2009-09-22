function [xTip,yTip,xCtr,yCtr] = get_head(im,head_length,threshold,visData,use_centroid)
% Finds coordinates for the rostrum and center of the cranium form an image
% of a back-lit fish.

if nargin < 4
    visData = 1;
    figure
end

%% Process the image

% Threshold & invert
bw = im2bw(im,threshold);
bw = ~bw;

% Find the areas of objects, pick the largest, remove others
L = bwlabel(bw,4);
stats = regionprops(L,'Area','Centroid');
maxArea = 1;
for i = 1:length(stats)
   maxArea = max([maxArea stats(i).Area]);
end
bw_head = bwareaopen(bw,maxArea,4);

% Process image to make the head a solid shape
se = strel('disk',2);
bw_head = imclose(bw_head,se);
bw_head = imfill(bw_head,'holes');

% Find head boundary, check for single image
[B,L2] = bwboundaries(bw_head,'noholes');
if length(B) > 1
    error('more than one body found in image');
else
    boundary = B{1};
end

% Remove coordinates at edges
idx = find(boundary(:,2)~=max(boundary(:,2)));
boundary = boundary(idx,:);
idx = find(boundary(:,2)~=min(boundary(:,2)));
boundary = boundary(idx,:);
idx = find(boundary(:,1)~=max(boundary(:,1)));
boundary = boundary(idx,:);
idx = find(boundary(:,1)~=min(boundary(:,1)));
boundary = boundary(idx,:);
clear idx

% Take only coordinates within a head length
%Note: assumes leftward directed swimming
maxX = min(boundary(:,2)) + head_length;
idx = find(boundary(:,2) < maxX);
boundary = boundary(idx,:);
x = boundary(:,2);
y = boundary(:,1);
clear idx

% Visualize steps in this transformation
if 0 % visData
    RGB = label2rgb(L);
    subplot(2,2,1)
    imshow(bw)
%     hold on
%     imshow(RGB)
%     hold off
    subplot(2,2,2)
    imshow(bw_head)
    hold on
    plot(x,y, 'r.')
    clear RGB
    hold off
end

clear i maxArea RGB idx B L L2 maxX stats
clear se threshold

%% Sort coordinates by radial position

    % Find center & radial position
    cntr = [mean(x) mean(y)];
    xT = x-cntr(1);
    yT = y-cntr(2);
    radT = atan2(yT,xT);
    
    % Find & sort positive radial coordinates
    radT_posI       = find(radT>0);  
    xTPos           = xT(radT_posI);
    yTPos           = yT(radT_posI);
    radT_pos        = radT(radT_posI);
    [radT_pos,indx] = sort(radT_pos);
    xTPos           = xTPos(indx)+cntr(1);
    yTPos           = yTPos(indx)+cntr(2);
    
    % Find & sort negative radial coordinates
    radT_negI       = find(radT<0);
    xTNeg           = xT(radT_negI);
    yTNeg           = yT(radT_negI);
    radT_neg        = radT(radT_negI);
    [radT_neg,indx] = sort(radT_neg);
    xTNeg           = xTNeg(indx)+cntr(1);
    yTNeg           = yTNeg(indx)+cntr(2);
    
    % Redefine x & y, based on assortment by radial position
    x = [xTPos; xTNeg];
    y = [yTPos; yTNeg];
    
    clear cntr xT yT radT radT_posI xTPos yTPos ratT_pos
    clear radT_pos indx xTPos yTPos radT_negI xTNeg yTNeg
    clear radT_neg radT_neg xTNeg yTNeg
    
%% Find center of head

if use_centroid
    leftofhead = [ones(size(bw_head,1),round(min(boundary(:,2))+head_length)) ...
        zeros(size(bw_head,1),...
        size(bw_head,2)-round(min(boundary(:,2))+head_length))];
    L = bwlabel(bw_head&leftofhead,4);
    stats = regionprops(L,'Centroid');
    
    xCtr = stats.Centroid(1);
    yCtr = stats.Centroid(2);
else
    xCtr = mean([x(1) x(end)]);
    yCtr = mean([y(1) y(end)]);
end

clear L stats boundary head_length

%% Calculate head curvature, select max curvature
    
    t = cumsum(sqrt([0,diff(x')].^2 + [0,diff(y')].^2))';
    ti = linspace(t(1),t(end),1000)';
    
    % Smoothing spline, calc coords, curvature (kappa)   
    tol = 1.e3;
    sp = spaps(t',[x y]',tol);
    dsp = fnder(sp); dspt = fnval(dsp,ti); 
    ddspt = fnval(fnder(dsp),ti);
    kappa = abs(dspt(1,:).*ddspt(2,:)-dspt(2,:).*ddspt(1,:))./...
          (sum(dspt.^2)).^(3/2);
    tmp     = fnval(sp,ti)';
    xS      = tmp(:,1);
    yS      = tmp(:,2);
    
    % Focus on the nose
    radT    = atan2(yS-yCtr,xS-xCtr);
    iNose   = find(abs(radT)>3*pi/4);
    kappa   = kappa(iNose);
    xS      = xS(iNose);
    yS      = yS(iNose);
    ti      = ti(iNose);
    
    %Find the coords for the max curvature around the nose
    iMax = find(kappa==max(kappa));
    xTip = xS(iMax);
    yTip = yS(iMax);
    
    if 1 %visData
%         subplot(2,2,3)
%         plot(x,y,'.',xS,yS,'r-')
%         hold on
%         plot([xCtr xTip],[yCtr yTip],'o-m')
%         axis equal
%         xlabel('x');ylabel('y')
        subplot(2,2,3)
        plot(ti,kappa)
        hold on
        plot([ti(iMax) ti(iMax)],[min(kappa) max(kappa)],'r-')
        xlabel('arclength'); ylabel('kappa')
        hold off
        subplot(2,2,4)
        imshow(im)
        hold on
        plot([xTip xCtr],[yTip yCtr],'r-',...
            xTip,yTip,'or')
        hold off
        pause(.01);
    end
    
    clear t ti xi yi tol sp dsp dspt ddspt kappa iMax tmp xS yS

    
