function calConstant = zCalibrateM(im)


f       = figure;
set(f,'Color','k');
set(f,'DoubleBuffer','on');

tVal = 0.5;n=0;i=1;x=[];y=[];
while 1 == 1
    imshow(im);hold on
    h = title(['Frame number = ' num2str(i)]);
    set(h,'Color','w');
    plot(x,y,'r+-')
    [xi,yi,but] = ginput(1);
    if isempty(but)
        break
    elseif but==30
        %increase frame
        if i == numFrames,
            i == numFrames;
        else
            i=i+1;
        end
        x = [];y = [];
    elseif but==46
        %increase frames by 10
        if i+10>numFrames
            i = numFrames;
        else 
            i = i+10;
        end
    elseif but==44
        %decrease frames by 10
        if i-10 <1
            i=1;
        else
            i = i-10;
        end
    elseif but==31
        %decrease frame
        if i==1
            i = 1;
        else
            i = i-1;
        end
        x = [];y = [];
    elseif but==102 | but==70
        %Choose first frame
        i = 1;
    elseif but==109 | but==77
        %Choose middle frame
        i = round(numFrames/2);
    elseif but==108 | but==76
        %Choose middle frame
        i = numFrames;
    elseif but==1
        n = n+1;
        x(n) = xi;
        y(n) = yi;
    elseif but==3
        n = n-1;
        x = x(1:n);
        y = y(1:n);
    end
end
close
if ~length(x)==2
    error(['Incorrect number of points!']);
end
dist        = (diff(x).^2 + diff(y).^2).^.5;
disp(' ');
knownDist   = inputdlg('What is the known distance (note: dont give units)?  ');

calConstant = str2num(knownDist{1}) ./ dist;