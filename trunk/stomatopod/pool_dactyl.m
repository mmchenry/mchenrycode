function pool_dactyl
% Pools data on the striking body from CT (simple_ct.m) and hydrodynamic 
% analyses (acquire.m)

% Density of water (kg m^-3)
rho = 1024;


%% Drag analysis

rPath = ['/Volumes/data_commuter/Projects/Patek_project/' ...
         'scaling data'];

% Get male data
a = dir([rPath '/males/smithii*']);

for i = 1:length(a)
    % Load drg
    load([rPath '/males/' a(i).name filesep 'drag_morph.mat'])
    
    % Calculate I for water
    I_w = 0.25.*pi.*rho.*trapz(drg.h(2:end),drg.r.^2 .* drg.T.^2);
    
    % Store
    males(i,1) = drg.D;
    I_water_male(i,1) = I_w;
    m_indiv(i,1) = str2num(a(i).name(end-2:end));
    
    clear I_w drg
end

% Display male results
[muhat,sigmahat,muci,sigmaci] = normfit(males);

disp(' ')
disp( 'D in males --------------- ')
disp(['   Mean = ' num2str(muhat)])
disp(['   L1   = ' num2str(muci(1))])
disp(['   L2   = ' num2str(muci(2))])
disp(' ')

%[m_indiv males]

clear muhat muci a


% Get female data
a = dir([rPath '/females/smithii*']);

for i = 1:length(a)
    % Load drg
    load([rPath '/females/' a(i).name filesep 'drag_morph.mat'])
    
     % Calculate I for water
    I_w = 0.25.*pi.*rho.*trapz(drg.h(2:end),drg.r.^2 .* drg.T.^2);
    
    % Store
    females(i,1) = drg.D;
    I_water_fe(i,1) = I_w;
    f_indiv(i,1) = str2num(a(i).name(end-2:end));
end

% Display female results
[muhat,sigmahat,muci,sigmaci] = normfit(females);

disp(' ')
disp( 'D in females --------------- ')
disp(['   Mean = ' num2str(muhat)])
disp(['   L1   = ' num2str(muci(1))])
disp(['   L2   = ' num2str(muci(2))])
disp(' ')

% Display pooled results
[muhat,sigmahat,muci,sigmaci] = normfit([males;females]);

disp(' ')
disp( 'D in both --------------- ')
disp(['   Mean = ' num2str(muhat)])
disp(['   L1   = ' num2str(muci(1))])
disp(['   L2   = ' num2str(muci(2))])
disp(' ')


% Pool and sort results
ID = [m_indiv;f_indiv];
Dall = [males;females];
Iall = [I_water_male;I_water_fe];

[ID,I] = sort(ID);

[ID Dall(I)]

Iall(I)

clear muhat muci a

clear all




%% CT Scan analysis

% Paths to the CT scans
dPath{1} = '/Volumes/data_commuter/Projects/Patek_project/ct_scans/16bitunsigned';
dPath{2} = '/Volumes/data_commuter/Projects/Patek_project/ct_scans/16bitunsigned052210smithiiUL1left';
dPath{3} = '/Volumes/data_commuter/Projects/Patek_project/ct_scans/16bitunsigned072910smithiiUL2left';

clrs{1} = 'r';
clrs{2} = 'b';
clrs{3} = 'g';

% Mass, in milligrams
mass(1,1) = 211.0;
mass(2,1) = 31.6;
mass(3,1) = 402.2;


figure

% Loop through dactyls
for i = 1:3
    
    % Load I data for each slice
    load([dPath{i} filesep 'I_by_slice'])
  
    if i==1 || i==3
        d.zVal = d.zVal(end:-1:1);
        
    end
    
    % Trim zeros pixel values from data
    idx = d.pixIntensity~=0;   

    d.zVal          = d.zVal(idx);
    d.pixIntensity  = d.pixIntensity(idx);
    d.pixr2         = d.pixr2(idx);
    
    % Define dimensionless position along dactyl
    h_rel = (d.zVal-min(d.zVal))./range(d.zVal);
    
    
    % Normalized I
    I_norm(i,1) = sum(d.pixr2) / sum(d.pixIntensity) / range(d.zVal)^2;
    r_max (i,1) = range(d.zVal)
    
    subplot(2,1,1)
    plot(d.zVal-min(d.zVal),d.pixIntensity,clrs{i})
    hold on
    
    subplot(2,1,2)
    plot(h_rel,d.pixr2./sum(d.pixIntensity)./range(d.zVal)^2,clrs{i})
    hold on
    
    clear d
end

% Add axis titles
subplot(2,1,2)
xlabel('Relative dactyl position')
ylabel('I*')

% Display pooled results
[muhat,sigmahat,muci,sigmaci] = normfit(I_norm);

disp(' ')
disp( 'I* --------------- ')
disp(['   Mean = ' num2str(muhat)])
disp(['   L1   = ' num2str(muci(1))])
disp(['   L2   = ' num2str(muci(2))])

I_norm
r_max.*1000

% Scaling of I/m  (should be pretty flat)
c = polyfit(log10(mass),log10(I_norm),1);

figure;
plot(log10(mass),log10(I_norm),'ob',...
     log10([mass(2) mass(3)]),log10([mass(2) mass(3)]).*c(1)+c(2),'-b')
title(['slope = ' num2str(c(1)) '   intercept = ' num2str(c(2))]);
axis square
xlabel('log10(mass)')
ylabel('log10(I/m)')