function p = get_params(indiv,p)
% Returns parameter values for requested individual

% Relative tolerence of the simulation
rel_tol = 1e-7;

% Code for each individual
indiv_num = [3 12 120 121 122 123 125 129 132 133 137];

% Linkage coordinates (m)
L1 = 10^-3 .* [5.95 6.07 6.536 7.598 6.533 6.885 7.104 7.276 7.798 7.023 5.75];
L2 = 10^-3 .* [3.755 3.903 4.227 4.871 4.1097 4.4355 4.324 4.812 4.9987 4.4675 3.627];
L3 = 10^-3 .* [0.425 0.4456 0.5676 0.595 0.6478 0.614 0.651 0.5832 0.7603 0.5871 0.4925];
L4 = 10^-3 .* [6.676 7.335 7.5588 8.596 7.511 7.988 7.91 8.293 8.718 7.839 6.467];

% % Length between points f and g (m)
% L5 = 10^-3 .* [6.058 7.149 6.958 7.994 6.784 6.964 7.494 7.254 7.527 ...
%                7.4586 5.973];
           
% Dactyl length (mm) -- calculated as the length btwn point b and the
% distal end of the dactyl
dac_len = 10^-3 .* [9.627 10.48 10.902 12.372 10.4549 11.179 11.357 ...
                    11.7348 12.105 11.439 9.4076];    
                
% Mass of the striking body (kg) -- estimated from scalign relationships of
% CT specimens
dac_mass =  10^-6 .* 10.^(3.13.*log10(dac_len.*1000) -0.885);              
           
% Moment of inertia for striking body (kg m^2) (calculated using I* = 0.28608)
dac_I = 0.28608 .* (( dac_len).^2) .* dac_mass;

% Length of COM of striking body
L5 = sqrt(dac_I./dac_mass);

% Distance between points A & F, determined by propodus thickness
h_AF  = 10^-3 .* [1.11 1.248 1.328 1.397 1.078 1.33 1.275 1.2856 ...
                  1.29 1.396 1.056];

% Initial angle of the dactyl (angle btwn link 4 and axis btwn f and g)
gamma = (pi/180).*[58.629 47.361 60.0725 57.274 56.697 57.166 54.821...
                   59.025 58.4864 56.037 58.978];

% Initial input angle (rad)
thetaStart = (pi/180).*[78.92 85.21 80 78.76 80 80.35 77.58 78.4 77.42 ...
                        76.63 78.39];

% Resting angle for torsion spring
thetaRest = (pi/180).*[92.18 99.44 92.31 91.76 92.49 92.92 93.41 90.63 ...
                       91.53 89.49 91.57];
                   
% Linear spring stiffness (N/m)               
k_lin = 10.^3 .*[43.777 41.376 27.322 33.649 23.271 72.302 38.244 23.122 ...
                 37.159 35.559 49.321];
             
% Angle between link 2 and wire for applying load during materials testing
lambda = (pi/180).*[26.868 27.258 25.787 28.78 28.125 25.847 29.678 ...
                    28.091 27.934 28.617 25.897];  

% Distance between load and point A:
dist_load = 10.^-3 .*[4.218 4.388 4.821 5.476 4.76 5.152 5.011 5.3895 ...
                      5.8346 5.0327 4.153];

% Torsion spring at mV joint (Nm/rad) 
k_spring = k_lin .* dist_load.^2 .* cos(lambda);

% Dimensionless drag index
D_indx = [0.0769 0.074 0.0756 0.0659 0.073 0.0621 0.0682 0.0741 0.0753 ...
          0.0743 0.0746];

% Moment of inertia for the water
waterI = 10^-8 .* [0.0335 0.0499 0.0539 0.0999 0.0556 0.0667 0.1002 ...
                   0.0981 0.1121 0.0753 0.0261];    
      
%TODO: Add values (and modify model) to deal with new drag coefficient

i = find(indiv==indiv_num);

if isempty(i)
    error('No matching individual')
end

p.L1 = L1(i);
p.L2 = L2(i);
p.L3 = L3(i);
p.L4 = L4(i);
p.L5 = L5(i);

p.h_AF  = h_AF(i);
p.gamma = gamma(i);

p.dacI    = dac_I(i);
p.dacMass = dac_mass(i);
p.dacLen  = dac_len(i);
p.D       = D_indx(i);

p.thetaStart = thetaStart(i);
p.thetaRest  = thetaRest(i);
p.kSpring = k_spring(i);

p.waterI = waterI(i);

p.rel_tol = rel_tol;

% Density of water (kg/m^3)
p.rho = 998;

