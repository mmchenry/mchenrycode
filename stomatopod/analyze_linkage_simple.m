function analyze_linkage
% Sensitivity analysis that examines KT and Range of Motion (ROM) for 
% different lengths in 4-bar linkage.


vis_pcolor = 1;

%% Prep to vary linkage dimensions

% Indivudal to simulate
indiv = 120;

kt_min = 5;
kt_max = 10;
rom_min = 10;
rom_max = 30;

% Load default geometry
p = get_params(indiv);

% Normalize by L1
L1_o = p.L1;
p.L1 = 1;
p.L2 = p.L2./L1_o;
p.L3 = p.L3./L1_o;
p.L4 = p.L4./L1_o;

numPts = 200;
lowerRange = .85;
upperRange = 1.15;

% 3 values to run each analysis (arranged in columns in final graph)
L1val     = [p.L1*lowerRange p.L1 p.L1*upperRange];
L2val     = [p.L2*lowerRange p.L2 p.L2*upperRange];
L3val     = [p.L3*lowerRange p.L3 p.L3*upperRange];
L4val     = [p.L4*lowerRange p.L4 p.L4*upperRange];


% Alter starting theta
%p.thetaStart = p.thetaStart.*lowerRange;
%p.thetaStart = p.thetaStart.*upperRange;

p_start = p;



%% Vary L2 at diff L3

L2M = linspace(L2val(1),L2val(3),numPts);

for i = 1:length(L3val)
    
    p.L3 = L3val(i);
    l_txt{i} = num2str(L3val(i));
    
    for j = 1:length(L2M)
        p.L2 = L2M(j);
        
        
        L = check_linkage(p,0,10^6);
        
        if isempty(L)
            KTvals(i,j) = nan;
            maxOutVals(i,j) = nan;
        else
            KTvals(i,j)     = max(L.KT_all);
            maxOutVals(i,j) = (L.thetaOutMax - L.thetaOutMin) * (180/pi);
        end
        
        
        clear L
    end
    
end


% Plot
figure;
subplot(2,2,1)
plot(L2M,KTvals)
xlabel('L2')
ylabel('max KT')
legend(l_txt,'Location','NorthWest')
axis square
ylim([4 10])
title('Varied L3')

subplot(2,2,3)
plot(L2M,maxOutVals)
xlabel('L2')
ylabel('range of motion')
axis square
ylim([12 28])


%% Vary L4 at diff L1

L4M = linspace(L4val(1),L4val(3),numPts);

for i = 1:length(L1val)
    
    p.L1 = L1val(i);
    l_txt{i} = num2str(L3val(i));
    
    for j = 1:length(L2M)
        p.L4 = L4M(j);
        
        L = check_linkage(p,0,10^6);
        
        if isempty(L)
            KTvals(i,j) = nan;
            maxOutVals(i,j) = nan;
        else
            KTvals(i,j)     = max(L.KT_all);
            maxOutVals(i,j) = (L.thetaOutMax - L.thetaOutMin) * (180/pi);
        end
        
        
        clear L
    end
    
end


% Plot

subplot(2,2,2)
plot(L4M,KTvals)
xlabel('L4')
ylabel('max KT')
legend(l_txt,'Location','NorthEast')
axis square
ylim([4 10])
title('Varied L1')

subplot(2,2,4)
plot(L4M,maxOutVals)
xlabel('L4')
ylabel('range of motion')
axis square
ylim([12 28])


