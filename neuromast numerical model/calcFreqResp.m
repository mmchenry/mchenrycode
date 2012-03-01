function c = calcFreqResp(c,h,relativeTo)
% Calculates the frequency response from the data in c for the bundle height.

if ~isfield(c,'M')
    error('need to run numerical solver first')
end

if nargin < 3
    relativeTo = 'hair bundle spd';
end

% grab vales from c
M       = c.M;
flwAmp = c.dispAmp;
freqs   = c.freqs;
heights = c.heights;


%Loop thru freqs, calculate performance at bundle
for ii = 1:size(M,2)
    omega       = 2*pi*freqs(ii);
    delta       = (2*c.mu/(omega*c.rho))^0.5;
    
    if strcmp(relativeTo,'hair bundle spd')
        stimFlw   = 1i*omega*flwAmp * (1-exp(-h*(1+1i)/delta));
        
    elseif strcmp(relativeTo,'hair bundle disp')
        stimFlw   = flwAmp * (1-exp(-h*(1+1i)/delta));
        
    elseif strcmp(relativeTo,'freestream disp')  
        stimFlw   = flwAmp;
        
    elseif strcmp(relativeTo,'freestream spd')    
        stimFlw   = 1i*omega*flwAmp;
        
    elseif strcmp(relativeTo,'nothing')       
        stimFlw   = 1;
        
    else
        error('Incorrect relativeTo input');
    end
    
    if c.numHeights==1
        bundleAmp = M(1,ii);
    else
        bundleAmp   = interp1(heights,M(:,ii),h);
    end
    bundleResp  = bundleAmp / stimFlw;
    
    
    c.sensitivity(ii,1) = abs(bundleResp);
    c.phase(ii,1)       = angle(bundleResp) / pi * 180;
    
    clear omega delta bundleSpd bundleAmp bundleResp
end



%calc sensitivity
c.peak_sense    = max(c.sensitivity);
c.peak_freq     = freqs(min(find(c.peak_sense==c.sensitivity)));