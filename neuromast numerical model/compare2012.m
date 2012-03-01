
function compare2012
% Evaluates different versions of the various models for the 2012 review
% with Sietse. 

%% Define what to run

comp_old = 0;

height_force = 0;

slide_vs_fixed = 1;


%% Define paths

% Define root path
root_path = '/Users/mmchenry/Dropbox/Chapter 5';

% Load 'm' structure (Sietse's data)
load([root_path filesep 'sietse_data.mat'])


%% Run numerical model

if comp_old

% Load Sieste's parameter values
c = sietse_params;

% Run simulation
c = numerical_twopart(c,'fixed','one part');
%c = numerical_twopart(c,'linear spring','one part');

% Calculated freq resp of sensitivity to freestream
c = calcFreqResp(c,c.bunHeight,'freestream spd');

end


%% Calc freq response from new method

if comp_old

% Calc flexural stiffness
EI = ((pi / 64) * c.baseDiameter^4) * c.E_matrix;

% Performa calculation individually for each freq
for i = 1:length(c.freqs)
    %v3d = interp1(c.heights,c.d3M(:,i),c.bunHeight);
    
    % Find third derivative of deflection close to base
    v3d = c.d3M(1,i);
    
    % Find force near base
    c_force = -EI*v3d;
    
    % Calculate flow speed (complex number)
    c_flowspd = 1i.*2.*pi.*c.freqs(i).*c.dispAmp;
    
    % Sensitivity (complex)
    sense_new(i,1) = c_force/(c.numHairs.*c.linStiff)/c_flowspd;
    
    % Clear variables for next loop
    clear c_force c_flowspd v3d
end

% Find amplitude and phase for sensitivity
c.amp_new = abs(sense_new);
c.ph_new  = angle(sense_new) / pi * 180;
    
clear EI sense_new

end


%% Plot Sietse's analytical prediction with the numerical model

if comp_old

figure

subplot(2,1,1)
loglog(c.freqs,c.sensitivity,'k',m.freq,m.amp_viaold,'r--')
legend('numerical','analytical')
xlabel('freq (Hz)')
ylabel('amp (s)')
title('Old method')

subplot(2,1,2)
semilogx(c.freqs,c.phase,'k',m.freq,m.ph_viaold,'r--')
xlabel('freq (Hz)')
ylabel('phase (deg)')

figure

subplot(2,1,1)
loglog(c.freqs,c.amp_new,'k',m.freq,m.amp_viaforce,'r--')
legend('numerical','analytical')
xlabel('freq (Hz)')
ylabel('amp (s)')
title('Force method')

subplot(2,1,2)
semilogx(c.freqs,c.ph_new,'k',m.freq,m.ph_viaforce,'r--')
xlabel('freq (Hz)')
ylabel('phase (deg)')

figure;

subplot(2,1,1)
loglog(c.freqs,c.sensitivity,'k',c.freqs,c.amp_new,'b--')
legend('old method','force method')
xlabel('freq (Hz)')
ylabel('amp (s)')
title('old vs. Force method')

subplot(2,1,2)
semilogx(c.freqs,c.phase,'k',c.freqs,c.ph_new,'b--')
xlabel('freq (Hz)')
ylabel('phase (deg)')

% Clear remnants of freq response analysis
clear c m

end


%% Calculate force as a function of height


if height_force

clrs = {'k' 'r' 'b'};

c = sietse_params;

c.freqs = [.03 .3 3]; % Hz

EI = ((pi / 64) * c.baseDiameter^4) * c.E_matrix;



% Run simulation
c = numerical_twopart(c,'fixed','one part');

c.force= -EI.*c.d3M;

%c.force= c.d3M;


figure
for i = 1:length(c.freqs)
    plot(abs(c.force(:,i)).*10^9,c.heights.*10^6,clrs{i})
    
    hold on
    
end

xlabel('Force (nN)')
%xlabel('3rd derivative of deflection ')
ylabel('height (microns)')
legend('.03 Hz','.3 Hz','3 Hz')

clear EI clrs i c

end


%% Run fixed and sliding base models

if slide_vs_fixed

% Load Sieste's parameter values
c = sietse_params;

% Run simulation
c_fixed = numerical_twopart(c,'fixed','one part');
c_slide = numerical_twopart(c,'linear spring','one part');

% Calculated freq resp of sensitivity to freestream
c_fixed = calcFreqResp(c_fixed,c.bunHeight,'freestream spd');
c_slide = calcFreqResp(c_slide,c.bunHeight,'freestream spd');

end


%% Force method of sensivity for both fixed and sliding models

if slide_vs_fixed
% Calc flexural stiffness
EI = ((pi / 64) * c.baseDiameter^4) * c.E_matrix;

% Performa calculation individually for each freq
for i = 1:length(c.freqs)
    %v3d = interp1(c.heights,c.d3M(:,i),c.bunHeight);
    
    % Find third derivative of deflection close to base
    v3d_fixed = c_fixed.d3M(1,i);
    v3d_slide = c_slide.d3M(1,i);
    
    % Find force near base
    c_force_fixed = -EI*v3d_fixed;
    c_force_slide = -EI*v3d_slide;
    
    % Calculate flow speed (complex number)
    c_flowspd = 1i.*2.*pi.*c.freqs(i).*c.dispAmp;
    
    % Sensitivity (complex)
    sense_fixed(i,1) = c_force_fixed/(c.numHairs.*c.linStiff)/c_flowspd;
    sense_slide(i,1) = c_force_slide/(c.numHairs.*c.linStiff)/c_flowspd;
    
    % Clear variables for next loop
    clear c_force c_flowspd v3d
end

% Find amplitude for sensitivity
c_fixed.amp_f = abs(sense_fixed);
c_slide.amp_f = abs(sense_slide);

% Find phase for sensitivity
c_fixed.ph_f  = angle(sense_fixed) / pi * 180;
c_slide.ph_f  = angle(sense_slide) / pi * 180;
    
clear EI sense_new 

figure;

subplot(2,2,1)
loglog(c.freqs,c_slide.sensitivity,'b',...
       c.freqs,c_slide.amp_f,'bo', ...
       c.freqs,c_fixed.sensitivity,'r--',...
       c.freqs,c_fixed.amp_f,'r.')
       
legend('slide, old','slide, force', 'fixed,old','fixed, force')
xlabel('freq (Hz)')
ylabel('amp (s)')
title('Normal bundle stiffness')

subplot(2,2,3)
semilogx(c.freqs,c_slide.phase,'b',...
         c.freqs,c_slide.ph_f,'bo',...
         c.freqs,c_fixed.phase,'r--',...
         c.freqs,c_fixed.ph_f,'r.')
xlabel('freq (Hz)')
ylabel('phase (deg)')

clear c

end


%% Run fixed and sliding base models (low linear stiffness)

if slide_vs_fixed

% Load Sieste's parameter values
c = sietse_params;

c.linStiff = c.linStiff./1000;

% Run simulation
c_fixed = numerical_twopart(c,'fixed','one part');
c_slide = numerical_twopart(c,'linear spring','one part');

% Calculated freq resp of sensitivity to freestream
c_fixed = calcFreqResp(c_fixed,c.bunHeight,'freestream spd');
c_slide = calcFreqResp(c_slide,c.bunHeight,'freestream spd');

end


%% Force method of sensivity for both (low linear stiffness)

if slide_vs_fixed
% Calc flexural stiffness
EI = ((pi / 64) * c.baseDiameter^4) * c.E_matrix;

% Performa calculation individually for each freq
for i = 1:length(c.freqs)
    %v3d = interp1(c.heights,c.d3M(:,i),c.bunHeight);
    
    % Find third derivative of deflection close to base
    v3d_fixed = c_fixed.d3M(1,i);
    v3d_slide = c_slide.d3M(1,i);
    
    % Find force near base
    c_force_fixed = -EI*v3d_fixed;
    c_force_slide = -EI*v3d_slide;
    
    % Calculate flow speed (complex number)
    c_flowspd = 1i.*2.*pi.*c.freqs(i).*c.dispAmp;
    
    % Sensitivity (complex)
    sense_fixed(i,1) = c_force_fixed/(c.numHairs.*c.linStiff)/c_flowspd;
    sense_slide(i,1) = c_force_slide/(c.numHairs.*c.linStiff)/c_flowspd;
    
    % Clear variables for next loop
    clear c_force c_flowspd v3d
end

% Find amplitude for sensitivity
c_fixed.amp_f = abs(sense_fixed);
c_slide.amp_f = abs(sense_slide);

% Find phase for sensitivity
c_fixed.ph_f  = angle(sense_fixed) / pi * 180;
c_slide.ph_f  = angle(sense_slide) / pi * 180;
    
clear EI sense_new 

subplot(2,2,2)
loglog(c.freqs,c_slide.sensitivity,'b',...
       c.freqs,c_slide.amp_f,'bo', ...
       c.freqs,c_fixed.sensitivity,'r--',...
       c.freqs,c_fixed.amp_f,'r.')
       
%legend('slide, old','slide, force', 'fixed,old','fixed, force')
xlabel('freq (Hz)')
ylabel('amp (s)')
title('Low bundle stiffness')

subplot(2,2,4)
semilogx(c.freqs,c_slide.phase,'b',...
         c.freqs,c_slide.ph_f,'bo',...
         c.freqs,c_fixed.phase,'r--',...
         c.freqs,c_fixed.ph_f,'r.')
xlabel('freq (Hz)')
ylabel('phase (deg)')

clear c

end



function c = sietse_params

% Parameters for all anlayses
 c.freqs         = [10.^linspace(-2,3,100)]';
 c.numHeights    = 50;  
 c.bunHeight     = 5.3e-6; %From Dinklo, 2005
 c.dispAmp       = 10 * 10^-6; %m
 c.E_matrix      = 80; %31 Pa
 c.EI_kino       = 2e-21; % 2e-21 N m^2
 c.bundleStiff   = 2.925e-14; %Nm/rad (van Netten & Kroese, 1987) 
 c.linStiff      = 0.001; %N/m (van Netten & Kroese, 1987) 
 c.rho           = 1000; %998 kg m^-3
 c.mu            = 0.001; %1e-3 Pa s

% Data from morphometric measurements (based on stiffness paper)
c.baseDiameter 	= 9e-6;
c.cupHeight     = 50e-6;
c.numHairs      = 10;

% The following should only matter if running the two-part model
c.midDiameter 	= 9e-6 ;
c.kinoHeight 	= 16e-6;



