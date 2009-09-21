function analyze_head_coords
% Runs analyze head code for coordinates hand selected from a video

%% Path definition

m_path = '/Volumes/Docs/Projects/head_swimming/erik_data';


%% Load coordinate data

load([m_path filesep 'coordinates'])


framerate = 500;
units = 'cm';

if ~(size(d.x,2)==3)
    error('The data need to be collected as 3 coordinates');
end

disp('NOTE: first point assumed to be the nose')

if max(isnan(d.x(:)))
    [iT,jT] = find(isnan(d.x),1,'first');
    disp(' ');
    warning(['Frame ' num2str(iT) ...
        ' has no entry -- using only frames prior to this one']);
    iRows   = 1:iT-1;
else
    iRows   = 1:size(d.x,1);
end

% xTip = d.x(iRows,1);
% yTip = d.y(iRows,1);
% xCtr = mean(d.x(iRows,2:3),2);
% yCtr = mean(d.y(iRows,2:3),2);
xTip = d.x(iRows,2);
yTip = d.y(iRows,2);
xCtr = mean(d.x(iRows,[1 3]),2);
yCtr = mean(d.y(iRows,[1 3]),2);
    
    
%% Transform subsequenct coordinates relative to first frame

% Define FoR from the head in the first frame
origin = [xCtr(1) yCtr(1)];
S      = localSystem(origin,[xTip(1) yTip(1)]);

% Execute transformation
ptsTip = globalToLocal([xTip yTip],origin,S);
ptsCtr = globalToLocal([xCtr yCtr],origin,S);

xT     = ptsTip(:,1);
yT     = ptsTip(:,2);
xC     = ptsCtr(:,1);
yC     = ptsCtr(:,2);

%clear ptsTip ptsCtr orgin S xTip xCtr yTip yCtr


%% Animate coordinates
if 1
    figure
    set(gcf,'DoubleBuffer','on')
    
    for i = 1:length(xT)
        plot([xCtr(i) xTip(i)],[yCtr(i) yTip(i)],'-',xCtr(i),yCtr(i),'o')
        %plot([xC(i) xT(i)],[yC(i) yT(i)],'-',xC(i),yC(i),'o')
        %axis([min(xC) max(xC) min(yC) max(yC)])
        axis([min(xCtr) max(xCtr) min(yCtr) max(yCtr)])
        title(['Frame ' num2str(i)])
        grid on
        pause(.1);
    end
   % close
end




%% FUNCTIONS ===================

function S = localSystem(P1,P2)
% Defines a transformation vector for a local coordinate system in an
% inertial frame of reference.  Uses P1 as the origin and P2 to find the
% direction of the x-axis.  Coordinates must be (1x2) vectors.

if size(P1,1)~=1 || size(P1,2)~=2 ||...
   size(P2,1)~=1 || size(P2,2)~=2
    error('Coordinates must be (1x2) vectors');
end

xAxis       = (P2-P1)./norm(P2-P1);
yAxis       = [xAxis(2);-xAxis(1)];
S           = [xAxis' yAxis];

function pts = localToGlobal(pts,origin,S)

if size(pts,2)~=2 || size(origin,2)~=2 
    error('Coordinates must be a (nx2) vector');
end

pts         = [inv(S)'*pts']';
pts(:,1)    = pts(:,1)+origin(1);
pts(:,2)    = pts(:,2)+origin(2);

function pts = globalToLocal(pts,origin,S)

if size(pts,2)~=2 || size(origin,2)~=2 
    error('Coordinates must be a (nx2) vector');
end

pts(:,1)    = pts(:,1)-origin(1);
pts(:,2)    = pts(:,2)-origin(2);
pts         = [S'*pts']';
