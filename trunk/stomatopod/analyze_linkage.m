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

numPts = 50;
lowerRange = .85;
upperRange = 1.15;

% 3 values to run each analysis (arranged in columns in final graph)
L2val     = [p.L2*lowerRange p.L2 p.L2*upperRange];
L3val     = [p.L3*lowerRange p.L3 p.L3*upperRange];
L4val     = [p.L4*lowerRange p.L4 p.L4*upperRange];


% Alter starting theta
%p.thetaStart = p.thetaStart.*lowerRange;
%p.thetaStart = p.thetaStart.*upperRange;

p_start = p;



%% Vary L2 & L3 at diff L4

[L2M,L3M] = giveValRangeMesh(p,'L2','L3',lowerRange,upperRange,numPts);

for k = 1:length(L4val)
    for i = 1:size(L2M,1)
        for j = 1:size(L2M,2)
            p.L2 = L2M(i,j);
            p.L3 = L3M(i,j);
            p.L4 = L4val(k);
            
            L = check_linkage(p,0);
            
            if isempty(L)
                KTvals(i,j,k) = nan;
                maxOutVals(i,j,k) = nan;
            else
                KTvals(i,j,k)     = max(L.KT_all);
                maxOutVals(i,j,k) = (L.thetaOutMax - L.thetaOutMin) * (180/pi);
            end
            clear L
        end
    end
end

% Plot
if vis_pcolor
    figure;
    for k = 1:length(L4val)
        subplot(2,length(L4val),k)
        h = pcolor(L2M,L3M,KTvals(:,:,k));
        set(h,'EdgeColor','none');
        set(h,'FaceColor','interp')
        axis square
        xlabel('L2')
        ylabel('L3')
        title('KT')
        colorbar('SouthOutside')
        caxis([kt_min kt_max])
        
        subplot(2,length(L4val),k+length(L4val))
        h = pcolor(L2M,L3M,maxOutVals(:,:,k));
        set(h,'EdgeColor','none');
        set(h,'FaceColor','interp')
        axis square
        xlabel('L2')
        ylabel('L3')
        title('ROM')
        colorbar('SouthOutside')
        caxis([rom_min rom_max])
    end
    
    subplot(2,3,1)
    title(['L4 = ' num2str(L4val(1)) '  KT'])
    subplot(2,3,2)
    title(['L4 = ' num2str(L4val(2)) '  KT'])
    subplot(2,3,3)
    title(['L4 = ' num2str(L4val(3)) '  KT'])
    
end

clear L2M L3M L4M KTvals maxOutVals



%% Vary L1 & L4 at diff L3

p = p_start;

[L1M,L4M] = giveValRangeMesh(p,'L1','L4',lowerRange,upperRange,numPts);

for k = 1:length(L3val)
    for i = 1:size(L1M,1)
        for j = 1:size(L1M,2)
            p.L1 = L1M(i,j);
            p.L4 = L4M(i,j);
            p.L3 = L3val(k);
            
            L = check_linkage(p,0);
            
            if isempty(L)
                KTvals(i,j,k) = nan;
                maxOutVals(i,j,k) = nan;
            else
                KTvals(i,j,k)     = max(L.KT_all);
                maxOutVals(i,j,k) = (L.thetaOutMax - L.thetaOutMin) * (180/pi);
            end
            clear L
        end
    end
end

if vis_pcolor
    figure;
    for k = 1:length(L3val)
        subplot(2,length(L3val),k)
        h = pcolor(L1M,L4M,KTvals(:,:,k));
        set(h,'EdgeColor','none');
        set(h,'FaceColor','interp')
        axis square
        xlabel('L1')
        ylabel('L4')
        title('KT')
        colorbar('SouthOutside')
        caxis([kt_min kt_max])
        
        subplot(2,length(L3val),k+length(L4val))
        h = pcolor(L1M,L4M,maxOutVals(:,:,k));
        set(h,'EdgeColor','none');
        set(h,'FaceColor','interp')
        axis square
        xlabel('L1')
        ylabel('L4')
        title('ROM')
        colorbar('SouthOutside')
        caxis([rom_min rom_max])
    end
    
    subplot(2,3,1)
    title(['L3 = ' num2str(L3val(1)) '  KT'])
    subplot(2,3,2)
    title(['L3 = ' num2str(L3val(2)) '  KT'])
    subplot(2,3,3)
    title(['L3 = ' num2str(L3val(3)) '  KT'])
    
end

clear L2M L3M L4M KTvals maxOutVals

return

%% Vary L3 & L4 at diff L2

p = p_start;

[L3M,L4M] = giveValRangeMesh(p,'L3','L4',lowerRange,upperRange,numPts);

for k = 1:length(L2val)
    for i = 1:size(L3M,1)
        for j = 1:size(L3M,2)
            p.L3 = L3M(i,j);
            p.L4 = L4M(i,j);
            p.L2 = L2val(k);
            
            L = check_linkage(p,0);
            
            if isempty(L)
                KTvals(i,j,k) = nan;
                maxOutVals(i,j,k) = nan;
            else
                KTvals(i,j,k)     = max(L.KT_all);
                maxOutVals(i,j,k) = (L.thetaOutMax - L.thetaOutMin) * (180/pi);
            end
            clear L
        end
    end
end

if vis_pcolor
    figure;
    for k = 1:length(L2val)
        subplot(2,length(L4val),k)
        h = pcolor(L3M,L4M,KTvals(:,:,k));
        set(h,'EdgeColor','none');
        set(h,'FaceColor','interp')
        axis square
        xlabel('L3')
        ylabel('L4')
        title('KT')
        colorbar('SouthOutside')
        caxis([kt_min kt_max])
        
        subplot(2,length(L4val),k+length(L4val))
        h = pcolor(L3M,L4M,maxOutVals(:,:,k));
        set(h,'EdgeColor','none');
        set(h,'FaceColor','interp')
        axis square
        xlabel('L3')
        ylabel('L4')
        title('ROM')
        colorbar('SouthOutside')
        caxis([rom_min rom_max])
    end
    
    subplot(2,3,1)
    title(['L2 = ' num2str(L2val(1)) '  KT'])
    subplot(2,3,2)
    title(['L2 = ' num2str(L2val(2)) '  KT'])
    subplot(2,3,3)
    title(['L2 = ' num2str(L2val(3)) '  KT'])
end

return


function [M1,M2] = giveValRangeMesh(p,field1,field2,lowVal,highVal,numVals)

if ~isfield(p,field1)
    error([field1 ' is not a field in p'])
end

if ~isfield(p,field2)
    error([field2 ' is not a field in p'])
end

lowEnd1  = lowVal .* eval(['p.' field1]);
highEnd1 = highVal .* eval(['p.' field1]);

vals1 = linspace(lowEnd1,highEnd1,numVals);

lowEnd2  = lowVal .* eval(['p.' field2]);
highEnd2 = highVal .* eval(['p.' field2]);

vals2 = linspace(lowEnd2,highEnd2,numVals);

[M1,M2] = meshgrid(vals1,vals2);

