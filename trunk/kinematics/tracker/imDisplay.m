function imBW = imDisplay(im,tVal)
%shows a bitmap image with a threshold level of pixels highlighted in color
imshow(im.cdata,im.colormap);
hold on;
fillColor           = [.43 .49 1];
imBW                = ~im2bw(im.cdata,im.colormap,tVal);
im.cdata(find(imBW))= 244.*ones(length((find(imBW))),1);
%cMap                = [[0:1./255:1]' [0:1./255:1]' [0:1./255:1]'];
cmap                = im.colormap;
cmap(245,:)         = fillColor;
imshow(im.cdata,cmap)
hold off;