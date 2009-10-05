function pivComparison(dataPath)
%Compares flow in the impulse tank between a piv measurement and the
%prediction from the law of continuity.
% Note: piv data and motor data assumed to be in the same directory.  Run
% "osiver" for the piv analysis.


% Default parameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ratio of working section area to piston area for continuity calculation
%areaRatio = 4.9; % Found from flow data
%areaRatio = 5.22; % Found from calipers
areaRatio = 1./(0.7595); % Found from calipers (new piston)

areaRatio = 1./(0.50);

% Prompt for data directory, if none given, load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    dataPath = uigetdir(pwd,'Choose directory that contains data files');
end
load([dataPath filesep 'piv_data.mat']);


% Run synchronizeData to align the timing of the data sources
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[amp,daq,vid]   = synchronizeData(dataPath);


% PIV data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([dataPath filesep 'piv_data.mat']);
[piv_t,piv_spd] = calcFlowSpeed(p);


% Sync PIV time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
piv_t =  vid.t(1) + piv_t;


% Flow speed predicted from continuity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
amp.U_pred = amp.vel ./ areaRatio; 


% Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
f = figure;
plot(piv_t.*1000-3,abs(piv_spd));
hold on
plot(amp.t_vel.*1000,amp.U_pred,'r')
legend('piv','prediction')
xlabel('time (ms)')
ylabel('flow speed (mm/s)')
grid on



% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [t,spd] = calcFlowSpeed(p)
%Finds the average speed across the working section for each interval

for i = 1:length(p)
    %Collect data for current frame
    grid_xsize = p(i).osiv_args.grid_xsize;
    
    x    = p(i).x;
    y    = p(i).y;
    u    = p(i).u;
    v    = p(i).v;
    t(i) = p(i).t;
    
    %Loop thru x positions in flow field
    for j = 1:grid_xsize
        xIndex  = j:grid_xsize:length(x)-(j-1);
        uMean(j)= mean(u(xIndex));
        vMean(j)= mean(v(xIndex));
        yMean(j)= mean(y(xIndex));
        xMean(j)= mean(x(xIndex));
        spd(j)  = (uMean(j)^2 + vMean(j)^2)^.5;
    end
    
    uSuper(i)  = mean(uMean);
    vSuper(i)  = mean(vMean);
    xSuper(i)  = mean(xMean);
    ySuper(i)  = mean(yMean);
    
    clear x y u v uMean vMean xMean yMean spd
end

spd      = (uSuper.^2 + vSuper.^2).^0.5;

