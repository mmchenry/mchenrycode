function [ra,rphi,ret,rv,rw,rz,rf] = nm4(f,z,E,r,h,k,n)
% function [rv,rw,ret,ra,rphi,rz] = nm3 (w,z,E,r,h)
% analytical single-beam superficial neuormast model
% Sietse van Netten, Primoz Pirih, 2012
%
% Slightly modified by Matt McHenry, 2012
%
%input parameters (independent variables)
% w   ... angular frequency of free-stream motion
% z   ... height of the calculation
% k   ... linear stiffness of single hair bundle
% n   ... number of hair bundles
%
% if w and z are not of the same size, then w*z combinations will be made
%
%parameters (scalars, to be fitted)
% E   ... cupular young module
% r   ... cupula radius
% h   ... cupula height
%
%returns
% rv ... complex angular velocity of the cupula at given (w,z)
% rw ... the angular frequency of free stream velocity
% ret .. a container with several calculation parameters
% ra ... real amplitude of the angular velocity (from rv)
% rphi . phase of the angular velocity (in degrees)
% rz ... the height of the cupula at the given (w,z)
%
% uses viscosity 1 mPoise, density 1kg/l, free stream velocity 1 m/s
% all units SI
%
%example
%[v w] = nm3; loglog(w,abs(v));

%these are the constant parameters
mu =  0.001; %0.001002;  %viscosity [Pa.s]
rho = 1000;	 %density [kg/m3]
u = 1; %free stream velocity; (m/s)

%add standard parameters if not given
if ~exist('E')
    E = 80; %Pa
end
if ~exist('r')
    r = 4.5e-6;
end
if ~exist('h')
    h = 50e-6;
end

%these are the standard testing ranges for our calculations
if ~exist('f')
    f = logspace(-1,4,501)';
end
if ~exist('z');
    z = linspace(1e-6,h,10);
    %z = 50e-6;
end
ret.z=z;

%convert from frequency to angular frequency
%match frequency and height input dimensions
%put those into the third dimension of the arrays

w=2.*pi.*f;
w=shiftdim(w); sw=size(w,1);
z=shiftdim(z); sz=size(z,1);
if sw~=sz
    w2=repmat(w,[1 sz]);
    z2=repmat(z',[sw 1]);
    sw2=size(w2);
    sz2=size(z2);
    w2=w2(:);
    z2=z2(:);
    w=shiftdim(w2);
    z=shiftdim(z2);
end
tw(1,1,:)=w; w=tw;
tz(1,1,:)=z; z=tz;

%boundary layer thickness
delta = sqrt (2*mu./rho./w); %dimension m;

%the dimensionless exponential
%Fz is the running variable used for water displacement at each height
%Fh is the cupular height used for d2 and d3 parameters
Fz = exp(-(1+1i).*z./delta);
Fh = exp(-(1+1i).*h./delta);

%the original from Matts paper
%calculate water velocity at the given z height
%vpar1 = u.*(1-exp(-z.*(1+i)./delta));
%we add -i/w to convert from velocity to displacement
%vpar1 = -i./w.*u.*(1-exp(-(1+1i).*z./delta));

%displacement of the fluid (dimensions m)
vpar1 = (-1i./w).*u.*(1-Fz);

% dimensionless
L = 0.5772 + log (r./sqrt(2)./delta);

% dimensionless
k = -L ./ (L.^2 + pi.*pi./16);

%bm the fluid damping coefficient unit (kg/m/s) or (Pa.s)
% 3 terms:
% (1) viscous,
% (2) inertial,
% (3) inertial of the added mass around cupula
bm = - 4.*pi.*mu.*k ...
    - 2i.*pi.*rho.*r.^2.*w ...
    + 1i.*pi.*pi.*mu.*k./L;

%I has dimensions of m^4
I = pi/4*r^4;

%wavenumber (fine)
kappa = (1i.*w.*bm./E./I).^(1/4); %dimension 1/m

%make factors a and b to be used for the matrix;
ta = exp(    kappa.*h); % dimensionless
tb = exp(1i.*kappa.*h); % dimensionless

%we need this for the matrix
%they are always in 3rd %dimension
a(1,1,:)=shiftdim(ta);
b(1,1,:)=shiftdim(tb);
one=ones(size(a));

%Moment vector:
%Used as the final dimension to build the displacement of the cupula
%Proportionality constants to build the final solution (also in m)
d2 = ( 2*u.*Fh        ) ./ (w .* delta.^2 .* kappa.^2); % dimension m
d3 = (-2*u.*(1+1i).*Fh) ./ (w .* delta.^3 .* kappa.^3); % dimension m
d4 = (  -u.*(1-1i)    ) ./ (w .* delta   .* kappa    ); % dimension m
%D has dimensions of m
D =  [zeros(size(d2)) ; d2 ; d3 ; d4];

%elements of M are dimensionless
M = [	one,  	one, 		one, 		one 		; ...
    a,      -b,		  1./a,	  -1./b	        ; ...
    a,      -1i.*b,	-1./a, 	1i./b		    ; ...
    one,	  +1i.*one, -one,	-1i.*one		];

if 1
    %SOLUTION FROM MAPLE
    Minv1_ = ...
        [[1i.*b.^2+2.*1i.*b.*a-b.^2+1+1i,            2.*1i.*b+1i.*a+b.^2.*a+1i.*b.^2.*a-a,     2.*1i.*b+a+1i.*b.^2.*a-b.^2.*a+1i.*a,   -1+1i+2.*1i.*b.*a+b.^2+1i.*b.^2           ]; ...
        [2.*1i.*b.*a+1i.*a.^2+a.^2-1+1i,              -1i.*a.^2.*b-2.*1i.*a+a.^2.*b-b-1i.*b,   -1i.*b-2.*a+1i.*a.^2.*b-a.^2.*b-b,      2.*a.*b+1-1i+a.^2+1i.*a.^2                ]; ...
        [(2.*1i.*b+1i.*a+b.^2.*a+1i.*b.^2.*a-a).*a,  (1i.*b.^2+2.*1i.*b.*a-b.^2+1+1i).*a,      -(-1+1i+2.*1i.*b.*a+b.^2+1i.*b.^2).*a,  -(2.*1i.*b+a+1i.*b.^2.*a-b.^2.*a+1i.*a).*a]; ...
        [(1i.*a.^2.*b+2.*1i.*a-a.^2.*b+b+1i.*b).*b,  -(2.*1i.*b.*a+1i.*a.^2+a.^2-1+1i).*b,     (2.*a.*b+1-1i+a.^2+1i.*a.^2).*b,        (-1i.*b-2.*a+1i.*a.^2.*b-a.^2.*b-b).*b   ]];
    %this is the divisor as copied from the Maple sheet;
    Minv2_ = -2/1i .* (b.^2+4.*a.*b + 1 + a.^2 + a.^2.*b.^2);
    %we must replicate it to match Minv1
    Minv2_ = repmat(Minv2_,[4,4,1]);
    Minv_ = Minv1_./Minv2_;
end

%comment: in the previous version i/2 was added at the end;
%Minv2 = 1 + a.^2 + b.^2 + 4.*a.*b + a.^2.*b.^2;
%Minv = -1i/2 .* Minv1 ./ Minv2;
if 0
    %SOLUTION FROM MAXIMA
    Minv1__ = [[(1i-1).*b.^2+2.*1i.*a.*b+1i+1, ...
        (1i+1).*a.*b.^2+2.*1i.*b+(1i-1).*a, ...
        (1i-1).*a.*b.^2+2.*1i.*b+(1i+1).*a, ...
        (1i+1).*b.^2+2.*1i.*a.*b+1i-1]; ...
        [2.*1i.*a.*b+(1i+1).*a.^2+1i-1, ...
        ((1-1i).*a.^2-1i-1).*b-2.*1i.*a, ...
        ((1i-1).*a.^2-1i-1).*b-2.*a, ...
        2.*a.*b+(1i+1).*a.^2-1i+1]; ...
        [(1i+1).*a.^2.*b.^2+2.*1i.*a.*b+(1i-1).*a.^2, ...
        (1i-1).*a.*b.^2+2.*1i.*a.^2.*b+(1i+1).*a, ...
        (-1i-1).*a.*b.^2-2.*1i.*a.^2.*b+(1-1i).*a, ...
        (1-1i).*a.^2.*b.^2-2.*1i.*a.*b+(-1i-1).*a.^2]; ...
        [((1i-1).*a.^2+1i+1).*b.^2+2.*1i.*a.*b, ...
        ((-1i-1).*a.^2-1i+1).*b-2.*1i.*a.*b.^2, ...
        2.*a.*b.^2+((1i+1).*a.^2-1i+1).*b, ...
        ((1i-1).*a.^2-1i-1).*b.^2-2.*a.*b]];
    %the divisor
    Minv2__ = (2.*1i.*a.^2+2.*1i).*b.^2+8.*1i.*a.*b+2.*1i.*a.^2+2.*1i;
    Minv2__ = repmat(Minv2__,[4,4,1]);
    Minv__ = Minv1__./Minv2__;
end

%Minv is dimensionless
Minv=Minv_;

% the four basic functions of the normal solution
% using the running variable for height
% W is dimensionless
W = [ exp(  1.* z .* kappa) ;
      exp( 1i.* z .* kappa) ;
      exp( -1.* z .* kappa) ;
      exp(-1i.* z .* kappa) ];

W_f = [ (kappa.^3 .* 1)  .* exp(  1.* z .* kappa) ;
        (kappa.^3 .* -1i).* exp( 1i.* z .* kappa) ;
        (kappa.^3 .* -1) .* exp( -1.* z .* kappa) ;
        (kappa.^3 .* 1i) .* exp(-1i.* z .* kappa) ];

%in the loop we calculate a hand-made inverse of the equation
%D = M * C ---->
%C = M\D;
%C = Minv * D;

%we go per each element in the third dimension (a combination of z and w)
for c=1:size(Minv,3);
    C(:,:,c) =     Minv(:,:,c) * D(:,:,c); % C has dimension meters,
    
    %vpar2(:,:,c) = 1;
    vpar2(:,:,c) = C(1,1,c) .* W(1,1,c) ...
        + C(2,1,c) .* W(2,1,c) ...
        + C(3,1,c) .* W(3,1,c) ...
        + C(4,1,c) .* W(4,1,c);
    
    % Calculate force on the cupula
%     v3par1 = (2*u*(1+1i))./(w(:,:,c).*delta(:,:,c).^3) .* ...
%                exp(-z(:,:,c).*(1+1i)./delta(:,:,c)) ;
%     
%     v3par2 =  C(1,1,c) .* W_f(1,1,c) ...
%             + C(2,1,c) .* W_f(2,1,c) ...
%             + C(3,1,c) .* W_f(3,1,c) ...
%             + C(4,1,c) .* W_f(4,1,c);
    
    %F(:,:,c)    =  E*I.* (v3par1 + v3par2);
    
    %v_bun(:,:,c) =  E*I.* (v3par1 + v3par2)./ (n*k);
    

    %numerical solution
    if 0
        C_(:,:,c) =    M(:,:,c)\D(:,:,c);
        vpar2_(:,:,c) =  C_(1,1,c) .* W(1,1,c) ...
            + C_(2,1,c) .* W(2,1,c) ...
            + C_(3,1,c) .* W(3,1,c) ...
            + C_(4,1,c) .* W(4,1,c);
    end
    
    %does not work: probably some fuckup with the inner product???
    %    vpar2__ (:,:,c)  = C(:,:,c)'    * W(:,:,c);
    %    vpar2___(:,:,c) = C_(:,:,c)'   * W(:,:,c);
end; %loop

% Remove nans (added by MJM)
vpar1(isnan(vpar1)) = 0;
vpar2(isnan(vpar2)) = 0;

%displacement v is a complex sum of water and cupular displacement
v    = vpar1 + vpar2 ;
%v_   = vpar1 + vpar2_ ;
%v__  = vpar1 + vpar2__ ;
%v___ = vpar1 + vpar2___ ;

ret.f=f;
ret.w=w;
%reshaping and returning values
rv=reshape(v,[sw sz]);
%rv_=reshape(v_,[sw sz]);
%rv__=reshape(v__,[sw sz]);
%rv___=reshape(v___,[sw sz]);

rf = 1; %reshape(v_bun,[sw sz]);

ra=abs(rv);
%ra_=abs(rv_);
%ra__=abs(rv__);
%ra___=abs(rv__);

rphi=angle(rv)*180/pi;
%rphi_=angle(rv_)*180/pi;
%rphi__=angle(rv__)*180/pi;
%rphi___=angle(rv___)*180/pi;

rz=reshape(z,[sw sz]);
rw=reshape(w,[sw sz]);

ret.E=E;
ret.r=r;
ret.h=h;

rp=sw2;

ret.v=reshape(v,rp);
%ret.v_=reshape(v_,rp);
ret.w=reshape(w,rp);
ret.vpar1=reshape(vpar1,rp);
ret.vpar2=reshape(vpar2,rp);
%ret.vpar2_=reshape(vpar2_,rp);

ret.ra = ra;
%ret.ra_ = ra_;

ret.rphi = rphi;
%ret.rphi_ = rphi_;

ret.delta=reshape(delta,rp);
ret.kappa=reshape(kappa,rp);

ret.Fz=reshape(Fz,rp);
ret.Fh=reshape(Fh,rp);
ret.W=reshape(W,[4 rp]);
ret.C=reshape(C,[4 rp]);
ret.D=reshape(D,[4 rp]);
ret.Minv=Minv;
ret.M=M;

%ret.v__=reshape(v__,rp);
%ret.v___=reshape(v___,rp);
%ret.vpar2__=reshape(vpar2_,rp);
%ret.vpar2___=reshape(vpar2_,rp);
%ret.rphi__ = rphi__;
%ret.rphi___ = rphi___;
%ret.ra__ = ra__;
%ret.ra___ = ra___;


end