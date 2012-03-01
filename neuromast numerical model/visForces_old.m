function visForces(c,freq)

%figPos = [1921         339        1200        1850];
figPos = [1         -16        1280         725];

if (min(freq) < min(c.freqs)) || (max(freq) > max(c.freqs))
    error(['requested freqeuncy out of range'])
end

% grab vales from c
bh      = c.bunHeight;
M       = c.M;
freqs   = c.freqs;
heights = c.heights;
L1      = c.kinoHeight;
dispAmp = c.dispAmp;
a1      = c.baseDiameter/2;
a2      = c.midDiameter/2;
indx1   = find(heights < L1);
indx2   = find(heights >= L1);
% Calculate 


% Step thru results to find solutions for requested freq values
for jj = 1:length(freq)
    currFreq        = freq(jj);
    omega           = 2 * pi * currFreq;
    
    % Boundary layer
    delta           = (2 * c.mu / c.rho / omega)^0.5;
    freeStrmSpd     = i * omega * dispAmp;   
    bLayerDisp(:,jj)    = dispAmp * (1 - exp(-heights * (1 + i)/delta));
    bLayerSpd(:,jj)     = i * omega .* bLayerDisp(:,jj);
    bLayerAccel(:,jj)   = -omega^2 .* bLayerDisp(:,jj);
 
    % Loop thru heights, evaluate at freq
    for ii = 1:length(heights)
        cupDisp(ii,jj)  = interp1(freqs,M(ii,:),currFreq);
        cupSpd(ii,jj)   = i*omega .* cupDisp(ii,jj);
        cupAccel(ii,jj) = -i*omega^2 .* cupDisp(ii,jj); 
    end
    
    % Fluid forces
    Re1     = omega * c.rho *a1^2 / c.mu;
    Re2     = omega * c.rho *a2^2 / c.mu;
    g1      = 0.577 + log(0.5*Re1^0.5);
    g2      = 0.577 + log(0.5*Re2^0.5);
    G1      = - g1 / (g1^2 + (pi/4)^2);
    G2      = - g2 / (g2^2 + (pi/4)^2);
    
    F_visc(:,jj)    = [4*pi*c.mu*G1 * ...
                       (bLayerSpd(indx1,jj) - cupSpd(indx1,jj));...
                       4*pi*c.mu*G2 * ...
                       (bLayerSpd(indx2,jj) - cupSpd(indx2,jj))];
    F_AR(:,jj)      = [(pi*c.rho*a1^2 - pi^2*c.mu*G1/g1/omega) .* ...
        		(bLayerAccel(indx1,jj) - cupAccel(indx1,jj)); ...
			(pi*c.rho*a2^2 - pi^2*c.mu*G2/g2/omega) .* ...
        		(bLayerAccel(indx2,jj) - cupAccel(indx2,jj))];
    F_buoy(:,jj)    = [c.rho*pi*a1^2 .* cupAccel(indx1,jj); ...
		       c.rho*pi*a2^2 .* cupAccel(indx2,jj)];

   % Structural forces
   E_matrix    = c.E_matrix;
   EI_kino     = c.EI_kino;
   I_base      = (pi / 64) * (2 * a1)^4;
   I_mid       = (pi / 64) * (2 * a2)^4;
   EI1         = (I_base * E_matrix) + (EI_kino * c.numHairs);
   EI2         = (I_mid * E_matrix);
   C_inert(:,jj)    = [pi*c.rho*a1^2 .* cupAccel(indx1,jj); ...
		       pi*c.rho*a2 .* cupAccel(indx2,jj)];
   %C_elast(:,jj)    = [EI1 .* 
end

% Plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
set(gcf,'Position',figPos);

%one row per frequency
for jj = 1:length(freq)
    numCols = 4;
    
    % Displacement plot
    subplot(length(freq),numCols,(jj-1)*numCols + 1)
    cAmp        = abs(cupDisp(:,jj));
    fSpd        = abs(bLayerSpd(:,jj));
    fSpd_norm   = fSpd ./ max(fSpd) .* max(cAmp);
    
    plot(cAmp,heights,'r-',fSpd_norm,heights,'b-')
    
    xlabel('Amplitude (m)')
    ylabel('height(m)')
    title(['Freq = ' num2str(freq(jj)) ' Hz'])
    if jj == 1, legend('cupula','norm flw');end
    
    clear cAmp fSpd fSpd_norm
    
    
    % Phase plot
    subplot(length(freq),numCols,(jj-1)*numCols + 2)
    cPh         = angle(cupDisp(:,jj)) ./ pi *180;
    blDisp      = angle(bLayerDisp(:,jj)) ./ pi *180;
    blSpd       = angle(bLayerSpd(:,jj)) ./ pi *180;
    blAccl      = angle(bLayerAccel(:,jj)) ./ pi *180;
    
    plot(cPh,heights,'r-');hold on
    plot(blDisp,heights,'b--')
    plot(blSpd,heights,'b-')
    plot(blAccl,heights,'b:')
    
    xlabel('Phase (deg)')
    if jj == 1, legend('cup disp','flw disp','flw spd','flw accel');end
    
    clear cPh blDisp blSpd blAccl
    
    
    % Force abs
    subplot(length(freq),numCols,(jj-1)*numCols + 3)
    Fv  = abs(F_visc(:,jj));
    Far = abs(F_AR(:,jj));
    Fb  = abs(F_buoy(:,jj));
    
    plot(Fv,heights,'r-'); hold on
    plot(Far,heights,'b-');
    plot(Fb,heights,'g-');
    
    xlabel('Force')
    if jj == 1, legend('F_visc','F_AR','F_buoy');end
    
    
    % Force phase
    subplot(length(freq),numCols,(jj-1)*numCols + 4)
    Fv  = 180 / pi .* angle(F_visc(:,jj));
    Far = 180 / pi .* angle(F_AR(:,jj));
    Fb  = 180 / pi .* angle(F_buoy(:,jj));
    
    plot(Fv,heights,'r-'); hold on
    plot(Far,heights,'b-');
    plot(Fb,heights,'g-');
    
    xlabel('Phase')
    if jj == 1, legend('F_visc','F_AR','F_buoy');end
    
end

