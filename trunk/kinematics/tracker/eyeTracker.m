function [d,xyData] = eyeTracker(mov,tVal,rEst,lEst,runInfo,roiSize,invert,...
                        subtract,startFrame,endFrame)
                    
roiType = 'fish eyes';

if runInfo.show
    figure;
    set(gcf,'Color','k'); 
    set(gcf,'DoubleBuffer','on');
    if runInfo.vidOut
        if mov.isTiff
            aviobj = avifile([mov.dirPath filesep mov.fileNames(1).name '_out'],'Compression','None','Quality',100,'fps',15);
        else
            aviobj = avifile([mov.dirPath filesep mov.fileName(1:end-4) '_out'],'Compression','None','Quality',100,'fps',15);
        end
    end
end

%Find intial coordinates and areas
im                  = grabFrame(mov,startFrame,invert,subtract);
[xNode,yNode]       = roiCoords(roiType,roiSize,lEst,rEst);
imROI               = roipoly(im.cdata,xNode,yNode);
imBW                = ~im2bw(im.cdata,im.colormap,tVal);
imBW                = imBW & imROI;
[rEst,lEst,rArea,lArea]... 
                    = giveTwoCoord(imBW,rEst,lEst); %get the likely coordinates for these eyes
%[coords,areas]      = giveSpotData(imBW);     
lAreaStart          = lArea;
lAreaOld            = lArea;
rAreaStart          = rArea;
rAreaOld            = rArea;

j                   = 1;
disp(' '); disp('started . . .');

%Loop through frames:
for i = startFrame:endFrame
    im                      = grabFrame(mov,i,invert,subtract);
    if j>2
        [rEst,lEst]         = giveEstimate(d.raw);
    end
    [xNode,yNode]           = roiCoords(roiType,roiSize,lEst,rEst);
    imROI                   = roipoly(im.cdata,xNode,yNode);
    imBW                    = ~im2bw(im.cdata,im.colormap,tVal);
    imBW                    = imBW & imROI;
    
    if ~max(imBW(:)), error('no spots are in the region of interest');end
        
    [imBW,tVal]             = giveTwoEyes(im,imROI,lAreaOld,rAreaOld,rAreaStart,lAreaStart,tVal,runInfo);
    
    %[rEst,rArea]        = giveTwoCoord(imBW,rEst); %get the likely coordinates for these eyes
    [rEst,lEst,rArea,lArea] = giveTwoCoord(imBW,rEst,lEst);
    [xNode,yNode]           = roiCoords(roiType,roiSize,lEst,rEst);
    
    d.raw.rEye.x(j,:)       = rEst(1);
    d.raw.rEye.y(j,:)       = rEst(2);
    d.raw.lEye.x(j,:)       = lEst(1);
    d.raw.lEye.y(j,:)       = lEst(2);
    
    xyData(j,:)             = [mean([rEst(1) lEst(1)]) mean([rEst(2) lEst(2)])];
        
    rAreaOld                = rArea;
    frameNum(j,:)           = i;
    j                       = j + 1;
    
    
  %Display:
    if runInfo.show
        %a                   = imDisplay(im,imBW,imROI,xNode,yNode,rEst,i,endFrame);
        imBW2                = imDisplay(im,tVal);  
        hold on
        plot(xNode,yNode,'y-',rEst(1),rEst(2),'ro',lEst(1),lEst(2),'sg')
        tit = title(['Frame ' num2str(i) ' out of ' num2str(endFrame)]); set(tit,'Color','w');
        hold off
        %disp(' ');
        %disp(['trackEyes:  FrameNum = ' num2str(i) ', FrameStop = ' num2str(endFrame)]);
        if runInfo.vidOut
            frmIm       = getframe(gca);
            aviobj      = addframe(aviobj,frmIm);
        end
        pause(.001)
    end
    
end


if runInfo.vidOut
    aviobj = close(aviobj);
end
if runInfo.show
    close
end
%save([mov.dirPath filesep 'xyData'],'d');
%csvwrite([mov.dirPath filesep 'xyData.txt'],d);
beep;disp('. . . done')


function [rEst,lEst] = giveEstimate(d)
cntrPnt         = [mean([d.lEye.x d.rEye.x],2) mean([d.lEye.y d.rEye.y],2)];
cntrPnt_disp    = [cntrPnt(end) - cntrPnt(end-1)];

ang2            = atan2([d.rEye.y(end) - cntrPnt(end,2)],[d.rEye.x(end) - cntrPnt(end,1)]);
ang1            = atan2([d.rEye.y(end-1) - cntrPnt(end-1,2)],[d.rEye.x(end-1) - cntrPnt(end-1,1)]);
ang_disp        = (ang2 - ang1);    

rEst            = [d.rEye.x(end) d.rEye.y(end)];
lEst            = [d.lEye.x(end) d.lEye.y(end)];

if (norm(cntrPnt_disp)>2) | ((180*ang_disp/pi)> 2)
    rEst_o    = rEst - cntrPnt(end,:);
    lEst_o    = lEst - cntrPnt(end,:);
    rEst(1)     =  [rEst_o(1) .* cos(ang_disp) + rEst_o(2) .* sin(ang_disp)] + cntrPnt(end,1);
    rEst(2)     =  [-rEst_o(1) .* sin(ang_disp) + rEst_o(2) .* cos(ang_disp)] + cntrPnt(end,2);
    lEst(1)     =  [lEst_o(1) .* cos(ang_disp) + lEst_o(2) .* sin(ang_disp)] + cntrPnt(end,1);
    lEst(2)     =  [-lEst_o(1) .* sin(ang_disp) + lEst_o(2) .* cos(ang_disp)] + cntrPnt(end,2);
end

function [imBW,lEst,rEst,area1,area2] = findEyes(im,imROI,tVal,lEst,rEst,lAreaOld,rAreaOld)
    imBW                = ~im2bw(im,tVal);
    imBW                = imBW & imROI;
    if  max(max(bwlabel(imBW)))==2
        imBW = imBW;
    else
        imBW = giveTwoEyes(im,imROI,lAreaOld,rAreaOld,tVal);
    end
    [rEst,lEst,rArea,lArea]...
                        = giveTwoCoord(imBW,rEst,lEst);

function [imBW,tVal] = giveTwoEyes(im,imROI,lAreaOld,rAreaOld,rAreaStart,lAreaStart,tVal,runInfo)
% Adjust tVal to keep area of eyes similar.  tValStep goes down an order of magnitude if 
% there is it misses the threshold in a step.
sDe         = 0;
sIn         = 0;
tValStep    = 10.^-2;
span        = 0.25;
while 1==1
    imBW                = ~im2bw(im.cdata,im.colormap,tVal);
    imBW                = imBW & imROI;
    [coords,areas]      = giveSpotData(imBW);
    [imBW,area1,area2]  = giveEyeImage(imBW,areas,lAreaOld,rAreaOld);
    if mean([area1 area2]) > (1+span) .*mean([lAreaStart rAreaStart])
        if sIn, tValStep = tValStep ./ 10; end
        tVal        = tVal - tValStep;
        if runInfo.show, disp(['      Decreasing threshold, tValStep =' num2str(tValStep)]);end
        sDe         = 1; 
        sIn         = 0;
    elseif mean([area1 area2]) < (1-span) .*mean([lAreaStart rAreaStart])
        if sDe, tValStep = tValStep ./ 10; end
        tVal        = tVal + tValStep;
        if runInfo.show, disp(['      Increasing threshold, tValStep =' num2str(tValStep)]);end
        sIn         = 1;
        sDe         = 0;
    else
        break
    end
    if log10(tValStep) < -5
        tValStep    = 10.^-2;
        span        = span + .05;
    end
end

function [imBW,area1,area2]  = giveEyeImage(imBW,areas,lAreaOld,rAreaOld);
imLabel         = bwlabel(imBW);
mArea           = mean([lAreaOld rAreaOld]);
areas           = [[1:length(areas)]' areas abs(areas-mArea)];
blobNum1        = min(areas(   find( min(areas(:,3))==areas(:,3) )   ,1));
area1           = areas( find(areas(:,1)==blobNum1) ,2);
im1             = imLabel == blobNum1;
areas           = areas(   find( ~(areas(:,1)==blobNum1) )  ,:);
blobNum2        = min(areas(   find( min(areas(:,3))==areas(:,3) )   ,1));
area2           = areas( find(areas(:,1)==blobNum2) ,2);
im2             = imLabel == blobNum2;
imBW            = im1 | im2;  

function [rEst,lEst,rArea,lArea] = giveTwoCoord(imBW,rEst,lEst)
imLabel         = bwlabel(imBW);
%First, grab centroids:
[coord1,area1]  = returnSpot(imLabel,imBW,1);
[coord2,area2]  = returnSpot(imLabel,imBW,2);
%Choose point most likely to correspond to each eye:
P1              = norm([coord1-rEst]) + norm([coord2-lEst]);
P2              = norm([coord2-rEst]) + norm([coord1-lEst]);
if P1 < P2
    rEst    = coord1;
    rArea   = area1;
    lEst    = coord2;
    lArea   = area2;
else   
    rEst    = coord2;
    rArea   = area2;
    lEst    = coord1;
    lArea   = area1;
end   

function [coords,areas]  = giveSpotData(imBW)
%Now, find the 2 most likely eye images, based on area:
imLabel     = bwlabel(imBW);
for i = 1:max(max(imLabel))
    [coordTemp,areaTemp]    = returnSpot(imLabel,imBW,i);
    coords(i,:)             = coordTemp;
    areas(i,:)              = areaTemp;
end

function [coord,areaS] = returnSpot(imLabel,imBW,lNum)
imSpot  = imBW & imLabel==lNum;
temp    = regionprops(bwlabel(imSpot),'Centroid','Area');
if isempty(temp)
    error('No spot in the region of interest');
end
coord   = temp.Centroid;
areaS   = temp.Area;

function a = imDisplay_old(im,imBW,imROI,xNode,yNode,rEye,frameNum,numFrames)
%imshow(im);hold on;
fillColor           = [.43 .49 1];
im.cdata(find(imBW&imROI))= 244.*ones(length((find(imBW&imROI))),1);
cMap                = [[0:1./255:1]' [0:1./255:1]' [0:1./255:1]'];
cMap(245,:)         = fillColor;
a                   = imshow(im.cdata,cMap);
hold on
plot(xNode,yNode,'g-');
eyeDisplay(rEye,frameNum,numFrames);
hold off

function eyeDisplay(rEye,frameNum,numFrames)
rSym = plot(rEye(1),rEye(2),'r+');
%set(rSym,'MarkerEdgeColor','k','MarkerFaceColor','k');
tit = title(['Frame ' num2str(frameNum) ' out of ' num2str(numFrames)]);
set(tit,'Color','w');

