function model_test
% Runs a series of simulations for individuals for which we have force data
% and compares prediction with measurement.  


% Create plot of individual simulations
visIndividual = 1;


%% Path definitions (computer specific)


% Path to text file of Mathematica commands
p.mathFile = ['"<</Volumes/data_commuter/Projects/Patek_project/m_files_podmodel/sim_code.txt"'];

% Path to the Mathematica kernel
p.kernelPath = '/Applications/Mathematica.app/Contents/MacOS/MathKernel';

% Root path for simulation data
p.simsPath = '/Volumes/data_commuter/Projects/Patek_project/sims';

% Path to force data
forcePath = '/Volumes/data_commuter/Projects/Patek_project/force_data';


%% Simulation parameters

% Precision of simulation (default = 10^-5, use 10^-7 for high precision)
p.maxError = 10^-7;

% Duration of simulation (s)
p.simDur    = 0.015;

% Time values to evaluate results(use 1000-5000)
p.t = linspace(0,p.simDur,1000);

% List of individuals for comparison
indivs = [120 121 122 123 125 129 132 133 137];


%% General parameter values
% Fixed for all individuals

% water I (kg m^2)
p.waterI      = 1e-8;

% 3rd moment of area of dactyl (m^5) (for 12 mm long dactyl, r = 3 mm)
p.dacA         = 0.25*(3e-3)*(12e-3)^4;

% Density of water (kg/m^3)
p.rho = 998;

%TODO: Adjust these to be custom for each individual


%% Load pooled data: fd

load([forcePath filesep 'force_pooled'])


%% Loop for each individual

for i = 1:length(indivs)
    
    % Current individual
    indiv = indivs(i);
    
    %% Extract measured data
    
    % Find data index for current indivdiual
    for j = 1:length(fd)
        
        if fd(j).indiv == indiv
            
            % Load t_start t_end
            load([forcePath filesep 'peak_' fd(j).filename '.mat'])
            
            % Calculate momentum
            p_cum = calc_momentum(fd(j).t,fd(j).F,pk_start,pk_end);
            
            r.p_meas(i)   = p_cum(end);
            r.maxForce(i) = max(fd(j).F);
            
            % Clear variables
            clear t_start t_end p_cum
            
            break
        end
        
        % Check for indivdiual match
        if j==length(fd)
            error('No match for individual in force data')
        end
    end

    
    %% Run simulations    
    
    % Report running
    
    % Individual-specific parameter values
    p = get_params(indiv,p);
    
    [d,result] = run_sim(p);
    
    r.p_model(i) = max(d.dacCOMSpd).*p.dacMass;
    r.dacI(i) = p.dacI;
    r.dacMass(i) = p.dacMass;
    
    clear d result
    
    disp(['Done ' num2str(i) ' of ' num2str(length(indivs))])
    
end

figure;
subplot(1,2,1)
boxplot(r.p_model)
ylim([0 7e-3])
title('model')
ylabel('Momentum')

subplot(1,2,2)
boxplot(r.p_meas)
ylim([0 7e-3])
title('measurement')


figure
subplot(3,1,1)
plot(r.p_model,r.p_meas,'o')
title('Momentum')
xlabel('Model')
ylabel('Measurement')
axis square
grid on

subplot(3,1,2)
plot(r.dacMass,r.maxForce,'o')
grid on 
axis square


return
if 0*visIndividual
    figure;
    subplot(3,1,1)
    plot(d.t.*1000,d.thetaIn.*180/pi,'r',d.t.*1000,d.thetaOut*180/pi,'b')
    legend('in ','out','Location','NorthWest');
    ylabel('angle (deg)')
    xlabel('time (ms)')
    grid on
    
    subplot(3,1,2)
    plot(d.t.*1000,(d.dacAngSpd.*180/pi)./1000)
    ylabel('SB angular spd (deg/ms)')
    xlabel('time (ms)')
    grid on
    
    subplot(3,1,3)
    E_tot = d.E_kin + d.E_drag + d.E_elastic;
    plot(d.t.*1000,d.E_elastic.*1000,'r',d.t.*1000,d.E_drag.*1000,'b',...
         d.t.*1000,d.E_kin.*1000,'g',d.t.*1000,E_tot.*1000,'k--')
    ylabel('Energy (mJ)')
    xlabel('time (ms)')
    grid on
    legend('elastic','drag','kinetic','total','Location','West');
    
    clear E_tot
end

%TODO: calc amd compare momentum values
%TODO: Set up to run all
%TODO: 




