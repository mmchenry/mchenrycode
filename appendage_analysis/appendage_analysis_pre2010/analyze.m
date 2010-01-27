fucntion analyze

% Run this after 'visualize'


% Calculate mass and morphology of segments along the smoothed body 
%--------------------------------------------------------------------------

% Check that the density vector has correct dimensions
if ~(size(morph.rho_tis)==size(w))
    error(['Tissue density vector must have the same dimensions '...
            'as the morphology']);
end
 
% Calculate area, volume, & mass of each element in the body
for i = 1:length(s)-1
    
    % Mean dimensions and density for segment
    ds                  = s(i+1)-s(1);
    a                   = mean([w(i) w(i+1)]) / 2;
    b                   = mean([h(i) h(i+1)]) / 2;
    rho                 = mean([morph.rho_tis(i) morph.rho_tis(i+1)]);
    
    % The circumerfence of an ellipse is approximated as . . .
    lambda              = (a-b)/(a+b);
    circ                = pi*(a+b)*(64-3*lambda^4)/(64-16*lambda^2);
    
    % Store results in morph
    morph.segXarea(i)   = 0.25*a*b;
    morph.segVol(i)     = 0.25*a*b*ds;
    morph.segMass(i)    = 0.25*a*b*ds*rho;
    morph.wetArea(i)    = circ*ds;
    
    clear ds a b rho lambda circ
end



% Calculate mechanical properties of the body
%--------------------------------------------------------------------------

% Calculate I


% Calculate the drag parameter


% Calculate the added mass parameter
