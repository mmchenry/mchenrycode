function [x,y] = roiCoords(varargin)
%Provides coordinates for a region of interest
action = varargin{1};


switch action
    case 'circle'
        %[x,y] = roiCoords('circle',roiDiameter,[xEstimate yEstimate]);
        roiRadius   = varargin{2}/2;
        rEst        = varargin{3};
        numPts      = 100;
        theta       = linspace(0,2*pi,numPts);
        x           = roiRadius .* cos(theta) + rEst(1);
        y           = roiRadius .* sin(theta) + rEst(2);
        
    case 'fish eyes'
        %[x,y] = roiCoords('fish eyes',roiWidth,[xR yR],[xL yL]);
        roiWidth    = varargin{2};
        lEst        = varargin{3};
        rEst        = varargin{4};
        zoneWidth   = 1 ;
        zoneHeight  = 1.25;
        zoneDepth   = 0.6;
        numPts      = 100;
        
        theta1      = [linspace(0,pi,round(numPts/2))]';
        theta2      = [linspace(pi,1.45*pi,round(numPts/2))]';
        theta3      = [linspace(1.55.*pi,2*pi,round(numPts/2))]';
        Xs          = [(zoneWidth/2) .* cos(theta1)];
        Ys          = [.9.*(zoneHeight/2) .* sin(theta1)];
        Xs          = [Xs; [(zoneWidth/2) .* cos(theta2)]];
        Ys          = [Ys; [(zoneDepth/2) .* sin(theta2)]];
        Xs          = [Xs; 0];
        Ys          = [Ys; (zoneDepth/2) .* sin(1.5.*pi).*0];%point
        Xs          = [Xs; [(zoneWidth/2) .* cos(theta3)]];
        Ys          = [Ys; [(zoneDepth/2) .* sin(theta3)]];
        
        origin      = [mean([rEst(1) lEst(1)]) mean([rEst(2) lEst(2)])];
        xAxis       = (rEst' - origin')./norm(rEst' - origin');
        yAxis       = [xAxis(2) -xAxis(1)]';
        S           = [xAxis yAxis];
        c           = [inv(S)'*[Xs Ys]']';
        x           = roiWidth .* c(:,1) + origin(1);
        y           = roiWidth .* c(:,2) + origin(2);
        
    case 'fish eyes with partition'
        %[x,y] = roiCoords('fish eyes',roiWidth,[xR yR],[xL yL]);
        roiWidth    = varargin{2};
        lEst        = varargin{3};
        rEst        = varargin{4};
        zoneWidth   = 1 ;
        zoneHeight  = 1.25;
        Xs          = [0    .1.*zoneWidth    .75.*zoneWidth    .75.*zoneWidth -.75.*zoneWidth -.75.*zoneWidth -.1.*zoneWidth 0]';
        Ys          = .75 .* [1.25.*zoneHeight   -.25.*zoneHeight   -.5.*zoneHeight...
            1.25.*zoneHeight 1.25.*zoneHeight -.5.*zoneHeight   -.25.*zoneHeight   1.25.*zoneHeight ]';
        %Ys          = [zoneHeight   -.3.*zoneHeight    -.3.*zoneHeight...
        %                1.25.*zoneHeight 1.25.*zoneHeight  -.3.*zoneHeight   -.3.*zoneHeight   .8.*zoneHeight ]';

        origin      = [mean([rEst(1) lEst(1)]) mean([rEst(2) lEst(2)])];
        xAxis       = (rEst' - origin')./norm(rEst' - origin');
        yAxis       = [xAxis(2) -xAxis(1)]';
        S           = [xAxis yAxis];
        c           = [inv(S)'*[Xs Ys]']';
        x           = roiWidth .* c(:,1) + origin(1);
        y           = roiWidth .* c(:,2) + origin(2);
end