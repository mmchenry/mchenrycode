function play


%% For kinematic analysis

%Use threshFinder to decide on a threshol value for this movie

mov_path = '/Volumes/Docs/Projects/head_swimming/s08_S0001';

head_length = 300;

fName = 'f1000001.tif';

im_path = [mov_path filesep fName];

im = imread(im_path);

roi.c = [86:975];
roi.r = [680:1209];

im = im(roi.r,roi.c);

threshold   = 0.44804;


get_head(im,head_length,threshold);








return



%% For morphological analysis

lat_path = '/Volumes/Docs/Projects/head_swimming/lat_head.jpg';
dors_path = '/Volumes/Docs/Projects/head_swimming/dors_head.jpg';

head_length = 400;


imL  = imread(lat_path);
imD = imread(dors_path);

[xD,yD] = head_profile(imD,head_length);
[xL,yL] = head_profile(imL,head_length);

cD = polyfit(yD,xD,2);
cL = polyfit(yL,xL,2);

figure;
imshow(imD);
hold on;
plot(polyval(cD,yD),yD,'r-')

figure;
imshow(imL);
hold on;
plot(polyval(cL,yL),yL,'r-')

