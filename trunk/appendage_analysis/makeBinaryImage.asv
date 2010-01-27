function makeBinaryImage(invertImage)
% Creates a binary image of a body of interest.
% Run this function before 'acquire' to create binary images to analyze.
% This works best on high contrast images 
% Set invertImage to 1 if the body you want ot analyze is lighter than its
% backdrop

if nargin < 1
    invertImage = 0;
end

if ~(invertImage==0 | invertImage==1)
    error('invertImage needs to be either 1 or 0');
end

% Browse to and load image
[fileName,filePath,indx]    = uigetfile('*.*','Choose an image file');
if indx==0  
    return
else
    im = imread([filePath filesep fileName]);
end

% Find first approximation for the threshold
g     = rgb2gray(im);
tVal  = graythresh(g);


% Interactively select a threshold value that works
f = figure;
subplot(2,1,1)
imshow(im)
subplot(2,1,2);
set(f,'DoubleBuffer','on');
bw    = chooseThresh(g,tVal,invertImage);
close

clear f tVal invertImage

% Cut blob of interest from connected blobs
bw = makeMask(bw,im);

clear im

% Select the blob that corresponds to the body of interest
bw = chooseBody(bw);
pause(.1)

% Save file
cd(filePath)
[fileName,savePath,indx] = uiputfile(fileName,'Save binary image');
if ~(indx==0)
    imwrite(bw,[savePath filesep fileName],'tiff','Compression','none');
end





%==========================FUNCTIONS=======================================


%--------------------------------------------------------------------------   
function bw = chooseThresh(img,tVal,invertImage)
%Used for finding coordinate points on a static image 'img'.

beep
disp(' ')
disp('Choose a threshold that fills in the shape of interest')
disp('  as a WHITE blob in a BLACK field, but do not allow');
disp('  the blob to fuse with any other blobs on the image')
disp(' ')
disp('  Up arrow -- increases the threshold')
disp('  Down arrow -- decreases the threshold')
disp('  = key -- micro increase in threshold')
disp('  - key -- micro decrease in threshold')
disp(' ')
disp('  Press return when done')
disp(' ')
    
% Loops each time a button is pushed, until return is pushed
but = 1;
while 1 == 1  
    % convert image to binary at threshold, invert image
    bw = im2bw(img,tVal);
    if ~invertImage
        bw      = ~bw;
    end
    
    % remove all objects containing fewer than 30 pixels
    bw = bwareaopen(bw,30);
    
    % fill gaps in white bodies
    se = strel('disk',2);
    bw = imclose(bw,se);
    
    % fill any large holes
    bw = imfill(bw,'holes');
    
    % display image
    imshow(bw);
    title('Make white blob on black field')
    xlabel(['tValue = ' num2str(tVal)])

    %prompt for input, adjust tVal accordingly
    [xi,yi,but] = ginput(1);
    if isempty(but)
        break
    elseif but == 30 %up arrow
        tVal = tVal + .01;
    elseif but == 31 %down arrow
        tVal = tVal - .01;
    elseif but == 61 %plus sign
        tVal = tVal + .001;
    elseif but == 45 %minus sign
        tVal = tVal - .001;    
    end
end 


%--------------------------------------------------------------------------   
function bw = chooseBody(bw)
% Returns binary image that includes only the blob chosen by the user

figure;
set(gcf,'DoubleBuffer','on');
imshow(bw)
title('Choose blob of interest');
[x,y,b] = ginput(1);
bw = bwselect(bw,x,y,4);
imshow(bw)
title('Final image')



%--------------------------------------------------------------------------   
function bw = makeMask(bw,im)
%Used for finding coordinate points on a static image 'bw'.

lWidth = 2;

subplot(2,1,2)
imshow(bw);
hold on;
subplot(2,1,1)
imshow(im)
title('Make polygon mask');
hold on

set(gcf,'DoubleBuffer','on');

disp(' '); disp(' ');
disp('Choose points to draw a mask on either image')
disp('Note: expand figure window for better precision')
disp(' ');
disp('Left mouse button picks points.');disp(' ');
disp('Right mouse button removes last point.');disp(' ');
disp('Press return when done.')

n = 0;
but = 1;
bwPoly = [];
while 1 == 1
    [xi,yi,but] = ginput(1);
    
    %Return pressed
    if isempty(but)
        break
        
    %Left click
    elseif but==1
        n = n+1;
        x(n) = xi;
        y(n) = yi;
             
        
    %Right click    
    elseif but==3
        if n-1 < 1
            n = 0;
            x = [];
            y = [];
        else
            n = n-1;
            x = x(1:n);
            y = y(1:n);
        end
    end   
    
    %Update plots
    subplot(2,1,2)
    bwPoly = roipoly(0,0,im,x,y);
    imshow(~bwPoly & bw)

    subplot(2,1,1)

    imshow(im)
    hold on
    fill(x,y,[1 0 0])
    hold off
    title('Make polygon mask');
end

if ~isempty(bwPoly)
    bw = ~bwPoly & bw;
end
close;
