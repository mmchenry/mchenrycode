function c = numerical_vibratingSphere(varargin)
% Find a numerical solution over m frequencies for n positions along the
% cupula.  Returns M, a m x n matrix of complex numbers describing the
% predicted deflection of the cupula.
% sFactor is a structure of scaling factor constants for length, time, and
% mass
%
% numerical_twopart(c) - runs model given c structure using default b
% conditions
%
% numerical_twopart(c,bConditions) - runs model given c structure using
% given b conditions
%
% The structure c describes the morphometrics of a cupula with the
% following fields:
% baseDiameter - diameter of the cupula at the base
% midDiameter - diameter of cupula at tips of kinocilia
% numHairs - number of hair cells (= number of kinocilia)
% kinoHeight - height of kinocilia
% cupHeight - total height of cupula
% freqs - the stimulus frequencies at which to evaluate 
% numHeights - the number of points along cupula to evaluate
% 
% bConditions - string specifying boundary condition type: |'torsion
% spring' |, 'fixed', or 'linear spring'
%
%

c = varargin{1};

if nargin > 1
    bConditions = varargin{2};
    if nargin > 2
        EIcondition = varargin{3};
    else
        EIcondition = 'two part';
    end
else
    bConditions = 'torsion spring';
    EIcondition = 'two part';
end


%Scaling constants (adjust to remove numerical instabilities)
% units.L         = 10^1; %in m  
% units.T         = 10^3; %in s  
% units.M         = 10^-3; %in kg 
units.L         = 10^-6; %in m  
units.T         = 10^-6; %in s  
units.M         = 10^-15; %in kg 


%Unit conversions
c.baseDiameter 	= c.baseDiameter  / units.L;
c.midDiameter 	= c.midDiameter  / units.L;
c.kinoHeight 	= c.kinoHeight / units.L;
c.cupHeight     = c.cupHeight / units.L;
c.freqs         = c.freqs .* units.T;
c.bunHeight     = c.bunHeight / units.L; 
c.dispAmp       = c.dispAmp / units.L; 
c.E_matrix      = c.E_matrix / units.M*units.L*units.T^2; 
c.EI_kino       = c.EI_kino / units.M/units.L^3*units.T^2; 
c.bundleStiff   = c.bundleStiff / units.M/units.L^2*units.T^2;  
c.linStiff      = c.linStiff / units.M *units.T^2;
c.rho           = c.rho / units.M*units.L^3; 
c.mu            = c.mu / units.M*units.L*units.T; 


%Define general parameters
E_matrix    = c.E_matrix;
EI_kino     = c.EI_kino;
bundleStiff = c.bundleStiff;
linStiff    = c.linStiff;
rho         = c.rho;
mu          = c.mu; 
dispAmp     = c.dispAmp;
L1          = c.kinoHeight;
L2          = c.cupHeight;
a1          = mean([c.baseDiameter c.midDiameter])/2;
a2          = c.midDiameter/2;
freqs       = c.freqs;

%Parameters for the sphere stimulus
spdSnd      = c.spdSnd;
p           = c.sphereDist; 
a           = c.sphereRadius;

I_base      = (pi / 64) * (2 * a1)^4;
I_mid       = (pi / 64) * (2 * a2)^4;

switch EIcondition
    case 'two part'
        EI1         = (I_base * E_matrix) + (EI_kino * c.numHairs);
        EI2         = (I_mid * E_matrix);
    case 'special'
        EI1         = (I_base * E_matrix) + (EI_kino * c.numHairs);
        EI2         = c.EIspecial / units.M/units.L^3*units.T^2;
        %EI2         = (I_mid * E_matrix);
    case 'one part'
        EI1         = (I_base * E_matrix) + (EI_kino * c.numHairs);
        EI2         = EI1;
end
clear EIcondition

L1 = c.kinoHeight;
L2 = c.cupHeight - c.kinoHeight;


torSpring   = bundleStiff .* c.numHairs;
linSpring   = linStiff .* c.numHairs;

if c.numHeights==1
    x1          = c.bunHeight;
    x2          = 0;
else
    num1        = round(c.numHeights .* c.kinoHeight / c.cupHeight);
    x1          = [linspace(0,c.kinoHeight, num1)]';
    x1          = x1(2:end);
    x2          = [linspace(c.kinoHeight,c.cupHeight,c.numHeights-num1+1)]';
    x2          = x2 - x2(1);
end

clear num1

%Define terms for governing equation that don't vary with freq
%  alph - (alpha) 2nd boundary condition - bundle stiffness
%  beta - (1 & 2) 5th & 6th boundary conditions - sets moments equal
%  phi - 7th boundary condition - sets moment at the tip
%  zta - (zeta) - 8th boundary condition - sets force at the tip
 
alph    = EI1 / torSpring;
beta1   = EI1;
beta2   = EI2;
gamma   = EI1 / linSpring;
phi     = 0;
zta     = 0; 

for ii = 1:length(freqs)
    
    % Define terms that vary with freq
    omega   = 2*pi*freqs(ii);
    Re1     = omega * rho *a1^2 / mu;
    Re2     = omega * rho *a2^2 / mu;
    g1      = 0.577 + log(0.5*Re1^0.5);
    g2      = 0.577 + log(0.5*Re2^0.5);
    G1      = - g1 / (g1^2 + (pi/4)^2);
    G2      = - g2 / (g2^2 + (pi/4)^2);
    
    %Stokes
    %b1      = omega * pi * ( - 4*mu*G1*i );
    %b2      = omega * pi * (- 4*mu*G2*i );
    b1      = omega * pi * (2*rho*a1^2*omega - 4*mu*G1*i - pi*mu*G1/g1);
    b2      = omega * pi * (2*rho*a2^2*omega - 4*mu*G2*i - pi*mu*G2/g2);
    
    
%     c.term1(ii) = 2*rho*a1^2*omega;
%     c.term2(ii) = - 4*mu*G2*i;
%     c.term3(ii) = - pi*mu*G2/g2;

    q_inf1  = - dispAmp * b1 / EI1; 
    q_inf2  = - dispAmp * b2 / EI2; 
    
    %max(abs(q_inf1(:)))
    
    k1      = (b1 / EI1)^0.25;
    k2      = (b2 / EI2)^0.25;
    delta   = ( 2 * mu / rho / omega )^0.5;
    
    s1      = -q_inf1/k1^4;
    s2      = -q_inf2/k2^4;
    cn      = -(1+i)/delta;

    t1      = (-s1/4) * ...
              ( i^0/(i^0-cn/k1) + ...
                i^1/(i^1-cn/k1) + ...
                i^2/(i^2-cn/k1) + ...
                i^3/(i^3-cn/k1) );
    t2      = (-s2/4) * (exp(cn*L1)) * ...
              ( i^0/(i^0-cn/k2) + ...
                i^1/(i^1-cn/k2) + ...
                i^2/(i^2-cn/k2) + ...
                i^3/(i^3-cn/k2) ); 
            
    % Find integration constants numerically:
    
    switch bConditions
        
        case 'torsion spring'
            coef    = find_coefs_torsionspring(k1,k2,L1,L2,alph,beta1,...
                        beta2,phi,zta,delta,s1,s2,cn,t1,t2);
        case 'linear spring'
            coef    = find_coefs_linearspring(k1,k2,L1,L2,alph,beta1,...
                        beta2,phi,zta,delta,s1,s2,cn,t1,t2,gamma);
        case 'fixed'
            coef    = find_coefs_fixed(k1,k2,L1,L2,alph,beta1,...
                        beta2,phi,zta,delta,s1,s2,cn,t1,t2);               
    end
    
    % Solution matrix:
    % Homogeneous part is the same as for a pressure field stimulus:
    
    fh1     = coef(1) * exp(k1 * x1) + ...
              coef(2) * exp(i * k1 * x1) + ...
              coef(3) * exp(-1 * k1 * x1) + ...
              coef(4) * exp(-i * k1 * x1);
    fh2     = coef(5) * exp(k2 * x2)  + ...
              coef(6) * exp(i * k2 * x2) + ...
              coef(7) * exp(-1 * k2 * x2) + ...
              coef(8) * exp(-i * k2 * x2); 
          
    % The nonhomogeneous part is unique to the flow stimulus, this version
    % works for a vibrating sphere:
    W       = dispAmp.*i.*omega;
    v       = mu/rho;
    fq1     = (i.*W)/(4.*omega) .* ...
                Rfunction(x1,p,k1,k2,spdSnd,v,omega,a,b1,EI1);
    fq2     = (i.*W)/(4.*omega) .* ...
                Rfunction(x2,p,k1,k2,spdSnd,v,omega,a,b2,EI2);
    %fq2     = s2 + t2 * exp(cn .* x2); 
    
    % General solution is the sum of the two parts:
    
    f1      = fh1 + fq1;
    f2      = fh2 + fq2;
    
    
    if c.numHeights==1
        c.M(1,ii)   = f1;
        c.dM(1,ii)  = Df1;
        c.d2M(1,ii) = DDf1;
        c.d3M(1,ii) = D3f1;
        c.d4M(1,ii) = D4f1;
    else
        c.M(:,ii)   = [f1(1:end-1); f2];
        c.dM(:,ii)  = [Df1(1:end-1); Df2];
        c.d2M(:,ii) = [DDf1(1:end-1); DDf2];
        c.d3M(:,ii) = [D3f1(1:end-1); D3f2];
        c.d4M(:,ii) = [D4f1(1:end-1); D4f2];
    end       
        
    clear fh1 fh2 fq1 fq2 f1 f2 Re1 Re2 omega
    
end

c.heights = [x1(1:end-1);x2+L1];



%Reverse unit conversions
c.baseDiameter 	= c.baseDiameter * units.L;
c.midDiameter 	= c.midDiameter * units.L;
c.kinoHeight 	= c.kinoHeight * units.L;
c.cupHeight     = c.cupHeight * units.L;
c.freqs         = c.freqs ./ units.T;
c.bunHeight     = c.bunHeight * units.L;
c.dispAmp       = c.dispAmp * units.L;
c.E_matrix      = c.E_matrix * units.M/units.L/units.T^2; 
c.EI_kino       = c.EI_kino * units.M*units.L^3/units.T^2; 
c.bundleStiff   = c.bundleStiff * units.M*units.L^2/units.T^2;  
c.linStiff      = c.linStiff * units.M /units.T^2;
c.rho           = c.rho * units.M/units.L^3; 
c.mu            = c.mu * units.M/units.L/units.T; 
c.M             = c.M .* units.L;
c.d2M           = c.d2M ./ units.L;
c.d3M           = c.d3M ./ units.L^2;
c.d4M           = c.d4M ./ units.L^3;
c.heights       = c.heights .* units.L;




function R = Rfunction(z,p,k1,k2,C,v,w,a,b,EI)
%This is flow for the profile normal to the direction of oscillation

% Define A1 and B1
% =========================================================
k = sqrt(w.^2 ./ (C.^2 + (4/3)*i*v.*w));
h = (1-i) * sqrt(w ./ (2*v));

A1 = - (3 + 3*i*h.*a - h.^2.*a.^2).*k.^3.*a.^3.*exp(i*k.*a) ./...
    (k.^2.*a.^2.*(1+i*h.*a) + (2+2*i*k.*a-k.^2.*a.^2).*h.^2.*a.^2);

B1 = (1/3) * (3 + 3*i*k.*a - k.^2.*a.^2).*h.^3.*a.^3.*exp(i*h.*a) ./...
    (k.^2.*a.^2.*(1+i*h.*a) + (2+2*i*k.*a-k.^2.*a.^2).*h.^2.*a.^2);

% Solve R
% =========================================================
beta = (b./EI).^(1/4);
R = -(i.*A1.*beta.^2)./(k1.^2) .* ...
            (lFunc(2,p-z,beta,k1)-lFunc(2,p+z,beta,k1))...
      ...                                                    
    - (A1.*beta.^3)./(k1.^3) .* ...
            (lFunc(3,p-z,beta,k1)-lFunc(3,p+z,beta,k1))...
      ...
    + (3.*B1.*beta)./k2 .* ...
            (lFunc(1,p-z,beta,k2)-lFunc(1,p+z,beta,k2))...
      ...
    - (3.*i.*B1.*beta.^2)./(k2.^2) .* ...
            (lFunc(2,p-z,beta,k2)-lFunc(2,p+z,beta,k2))...
      ...
    - (3.*B1.*beta.^3)./(k2.^3) .* ...
            (lFunc(3,p-z,beta,k2)-lFunc(3,p+z,beta,k2));
        

function lambda = lFunc(L,z,beta,alpha)
terms = 0;
for n = 0:3
    terms = terms + ...
        (-1)^(L+1) * i^(n*L) .* exp(i^n .* beta .* z) .* ...
        expint((i^n * beta + alpha) .* z)
end
lambda = (1./factorial(L-1)) .* terms;
            

