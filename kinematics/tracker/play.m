
%PLAY Summary of this function goes here
%   Detailed explanation goes here


%dev_info = imaqhwinfo('winvideo',1)
vid = videoinput('winvideo',1,'Y800_640x480');
%preview(vid)
% Set video input object properties for this application.
% Note that example uses both SET method and dot notation method.
set(vid,'TriggerRepeat',Inf);
vid.FrameGrabInterval = 5;

% Set value of a video source object property.
vid_src = getselectedsource(vid);
set(vid_src,'Tag','motion detection setup');

% Create a figure window.
%figure; 

% Start acquiring frames.
start(vid)

i = 1;
while(vid.FramesAcquired<=10) % Stop after 10 frames
    i               = i+1;
    [data,time]     = getdata(vid,1); 
    frameData(:,:,i)   = data(:,:,1,1);
    times(i)        = time;
    %im=data(:,:,:,1);
    %diff_im = imabsdiff(data(:,:,:,1),data(:,:,:,2));
    %image(im)
end
stop(vid)
delete(vid)

figure;
for i=1:size(frameData,3)
    imshow(frameData(:,:,i));
    pause(.5);
end
