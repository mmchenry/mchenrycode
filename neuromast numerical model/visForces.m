function visForces(c,freq)

%figPos = [1921         339        1200        1850];
%figPos = [1         -16        1280         725];

if (min(freq) < min(c.freqs)) || (max(freq) > max(c.freqs))
    error(['requested freqeuncy out of range'])
end

% grab vales from c
bh      = c.bunHeight;
M       = c.M;
D4M     = c.d4M;
freqs   = c.freqs;
heights = c.heights;
L1      = c.kinoHeight;
dispAmp = c.dispAmp;
a1      = mean([c.baseDiameter c.midDiameter])/2;
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
        cupAccel(ii,jj) = -omega^2 .* cupDisp(ii,jj); 
        cupD4M(ii,jj)   = interp1(freqs,D4M(ii,:),currFreq);
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
                    
    F_buoy(:,jj)    = [c.rho*pi*a1^2 .* bLayerAccel(indx1,jj); ...
                        c.rho*pi*a2^2 .* bLayerAccel(indx2,jj)];

   F_viscPH(:,jj)   = angle(F_visc(:,jj)/freeStrmSpd) * 180/pi;
   F_ARPH(:,jj)     = angle(F_AR(:,jj)/freeStrmSpd) * 180/pi;
   F_buoyPH(:,jj)   = angle(F_buoy(:,jj)/freeStrmSpd) * 180/pi;
   
   
   % Structural forces
   E_matrix    = c.E_matrix;
   EI_kino     = c.EI_kino;
   I_base      = (pi / 64) * (2 * a1)^4;
   I_mid       = (pi / 64) * (2 * a2)^4;
   EI1         = (I_base * E_matrix) + (EI_kino * c.numHairs);
   EI2         = (I_mid * E_matrix);
   
   C_inert(:,jj)    = [pi*c.rho*a1^2 .* cupAccel(indx1,jj); ...
                        pi*c.rho*a2^2 .* cupAccel(indx2,jj)];
   %C_elast(:,jj)    = F_visc(:,jj) + F_AR(:,jj) + F_buoy(:,jj) - C_inert(:,jj);
   C_elast(:,jj)      = [EI1 .*cupD4M(indx1,jj); EI2 .*cupD4M(indx2,jj)];
   C_inertPH(:,jj)    = angle(C_inert(:,jj)/freeStrmSpd) * 180/pi;
   C_elastPH(:,jj)    = angle(C_elast(:,jj)/freeStrmSpd) * 180/pi;
end

% Plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
%set(gcf,'Position',figPos);

%DEFLECTION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%one row per frequency
for jj = 1:length(freq)
    currFreq        = freq(jj);
    omega           = 2 * pi * currFreq;
    freeStrmSpd     = i * omega * dispAmp; 
    
    numCols = 2;
    
    % Displacement plot
    subplot(length(freq),numCols,(jj-1)*numCols + 1)
    cAmp        = abs(cupDisp(:,jj));
    fSpd        = abs(bLayerSpd(:,jj));
    fSpd_norm   = fSpd ./ max(fSpd) .* max(cAmp);
    
    plot(cAmp,heights,'r-',fSpd_norm,heights,'b-')
    axis square
    set(gca,'TickDir','out')
    set(gca,'TickLength',[.02 .02])
    set(gca,'YLim',[0 max(heights)])
    set(gca,'YTick',[0:1e-5:max(heights)])
    tmp=get(gca,'XLim');
    set(gca,'XLim',[tmp(1)-range(tmp)/20 tmp(2)]); clear tmp
    
    xlabel('Amplitude (m)')
    ylabel('height(m)')
    title(['Freq = ' num2str(freq(jj)) ' Hz'])
    if jj == 1, legend('cupula','norm flw');end
    
    clear cAmp fSpd fSpd_norm
    
    
    % Phase plot
    subplot(length(freq),numCols,(jj-1)*numCols + 2)
    cPh         = angle(cupDisp(:,jj)/freeStrmSpd) ./ pi *180;
    blDisp      = angle(bLayerDisp(:,jj)/freeStrmSpd) ./ pi *180;
    blSpd       = angle(bLayerSpd(:,jj)/freeStrmSpd) ./ pi *180;
    blAccl      = angle(bLayerAccel(:,jj)/freeStrmSpd) ./ pi *180;
    
    plot(cPh,heights,'r-');hold on
    plot(blDisp,heights,'b--')
    plot(blSpd,heights,'b-')
    %plot(blAccl,heights,'b:')
    axis square
    set(gca,'TickDir','out')
    set(gca,'TickLength',[.02 .02])
    set(gca,'XTick',[-90:45:90])
    set(gca,'XLim',[-90 90])
    set(gca,'YLim',[0 max(heights)])
    set(gca,'YTick',[0:1e-5:max(heights)])
    tmp=get(gca,'XLim');
    set(gca,'XLim',[tmp(1)-range(tmp)/20 tmp(2)+range(tmp)/20]); clear tmp
    
    xlabel('Phase (deg)')
    if jj == 1, legend('cup disp','flw disp','flw spd','flw accel');end
    
    clear cPh blDisp blSpd blAccl
    
    
end


%STRUCTURAL FORCES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for jj = 1:length(freq)
    numCols = 2;
    
    % Force abs
    subplot(length(freq),numCols,(jj-1)*numCols + numCols-1)
    C_e  = abs(C_elast(:,jj));
    C_i = abs(C_inert(:,jj));
    
    plot(C_e,heights,'r-'); hold on
    plot(C_i,heights,'r--');
    set(gca,'TickDir','out')
    set(gca,'TickLength',[.02 .02])
    set(gca,'YLim',[0 max(heights)])
    set(gca,'YTick',[0:1e-5:max(heights)])
    tmp=get(gca,'XLim');
    set(gca,'XLim',[tmp(1)-range(tmp)/20 tmp(2)]); clear tmp
    axis square
    xlabel('Force')
    if jj == 1, legend('C_e','C_i');end
    
    
%     % Force phase
%     subplot(length(freq),numCols,(jj-1)*numCols + numCols)
%     C_e  = abs(C_elastPH(:,jj));
%     C_i = abs(C_inertPH(:,jj));
%     
%     plot(C_e,heights,'r-'); hold on
%     plot(C_i,heights,'b-');
%     set(gca,'TickDir','out')
%     set(gca,'TickLength',[.02 .02])
%     set(gca,'YLim',[0 max(heights)])
%     set(gca,'YTick',[0:1e-5:max(heights)])
%     tmp=get(gca,'XLim');
%     set(gca,'XLim',[tmp(1)-range(tmp)/20 tmp(2)]); clear tmp
%     axis square
%     xlabel('Phase')
    
end


%FLUID FORCES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure;
%one row per frequency
for jj = 1:length(freq)
    numCols = 2;
    
    % Force abs
    subplot(length(freq),numCols,(jj-1)*numCols + numCols-1)
    Fv  = abs(F_visc(:,jj));
    Far = abs(F_AR(:,jj));
    Fbuoy = abs(F_buoy(:,jj));
    
    plot(Fv,heights,'b--'); hold on
    plot(Far,heights,'b-');
    plot(Fbuoy,heights,'k-');
    set(gca,'TickDir','out')
    set(gca,'TickLength',[.02 .02])
    set(gca,'YLim',[0 max(heights)])
    set(gca,'YTick',[0:1e-5:max(heights)])
    tmp=get(gca,'XLim');
    set(gca,'XLim',[tmp(1)-range(tmp)/20 tmp(2)]); clear tmp
    axis square
    xlabel('Force')
    if jj == 1, legend('F_v_i_s_c','F_A_R');end
    
    
    % Force phase
%     subplot(length(freq),numCols,(jj-1)*numCols + numCols)
%     Fv  = F_viscPH(:,jj);
%     Far = F_ARPH(:,jj);
%     
%     plot(Fv,heights,'r-'); hold on
%     plot(Far,heights,'b-');
%     set(gca,'TickDir','out')
%     set(gca,'TickLength',[.02 .02])
%     set(gca,'XTick',[-180:90:180])
%     set(gca,'XLim',[-180 180])
%     set(gca,'YLim',[0 max(heights)])
%     set(gca,'YTick',[0:1e-5:max(heights)])
%     tmp=get(gca,'XLim');
%     set(gca,'XLim',[tmp(1)-range(tmp)/20 tmp(2)+range(tmp)/20]); clear tmp
%     axis square
%     xlabel('Phase')
%     if jj == 1, legend('F_v_i_s_c','F_A_R');end
    
end


return

%STRUCTURAL & FLUID FORCES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for jj = 1:length(freq)
    numCols = 1;
    
    % Force abs
    subplot(length(freq),1,jj)
     C_tot = real(C_inert(:,jj)) + real(C_elast(:,jj));
     F_tot = real(F_visc(:,jj)) + real(F_AR(:,jj));
    
    %C_tot = real(C_elast(:,jj));
    %F_tot = real(F_visc(:,jj)) + real(F_AR(:,jj)) ;
    
    plot(C_tot,heights,'r-'); hold on
    plot(F_tot,heights,'b--');
    set(gca,'TickDir','out')
    set(gca,'TickLength',[.02 .02])
    set(gca,'YLim',[0 max(heights)])
    set(gca,'YTick',[0:1e-5:max(heights)])
    tmp=get(gca,'XLim');
    set(gca,'XLim',[tmp(1)-range(tmp)/20 tmp(2)]); clear tmp
    axis square
    xlabel('Force')
    if jj == 1, legend('C_t_o_t','F_t_o_t');end
    
   
end