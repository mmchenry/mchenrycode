function bw = changePixels(bw,x,y)
% Turns pixels black in binary at overlapping x and y values

offst = 0;
if length(x)>1
    for j = 1:2
        xPair = [x(end-1) x(end)] + offst;
        yPair = [y(end-1) y(end)];

        if xPair(1)==xPair(2)
            xPair(2) = xPair(2) + 3;
        end
        
        iMin = min(find(xPair==min(xPair)));
        if iMin==1, iMax=2;else, iMax=1;end
        xVals = ceil(xPair(iMin)):floor(xPair(iMax));
        
        yVals = round(interp1(xPair,yPair,xVals));

        % Turn overlapping pixels black
        for i = 1:length(xVals)
            bw(yVals(i),xVals(i))=0;
        end
        
        offst = offst + 1;
    end
else
    bw = bw;
end