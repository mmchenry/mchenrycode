function r = find_coefs(k1,k2,L1,L2,alph,beta1,beta2,phi,...
                            zta,delta,s1,s2,c,t1,t2)
% Finds numerical values for the integration constants of the cupula model
% with two parts and a torsion spring at the base.

Z(2,1)      = k1 - alph * k1^2;
Z(2,2)      = k1*i + alph * k1^2;
Z(2,3)      = -k1 - alph*k1^2;
Z(2,4)      = -i*k1 + k1^2*alph;
%Z(2,:) = [1-k1*alph, i+k1*alph, -1-k1*alph, -i+k1*alph]
Z(2,5:8)    = [0 0 0 0];

Z(1,:)      = [1 1 1 1 0 0 0 0];

Z(3,1)      = exp(k1*L1);
Z(3,2)      = exp(i*k1*L1);
Z(3,3)      = exp(-k1*L1);
Z(3,4)      = exp(-i*k1*L1);
Z(3,5:8)    = [-1 -1 -1 -1];

Z(4,1)      = k1*exp(k1*L1);
Z(4,2)      = i*k1*exp(i*k1*L1);
Z(4,3)      = -k1*exp(-k1*L1);
Z(4,4)      = -i*k1*exp(-i*k1*L1);
Z(4,5)      = -k2;
Z(4,6)      = -i*k2;
Z(4,7)      = k2;
Z(4,8)      = i*k2;

Z(5,1)      = beta1*k1^2*exp(k1*L1);
Z(5,2)      = -beta1*k1^2*exp(i*k1*L1);
Z(5,3)      = beta1*k1^2*exp(-k1*L1);
Z(5,4)      = -beta1*k1^2*exp(-i*k1*L1);
Z(5,5)      = -beta2*k2^2;
Z(5,6)      = beta2*k2^2;
Z(5,7)      = -beta2*k2^2;
Z(5,8)      = beta2*k2^2;

Z(6,1)      = beta1*k1^3*exp(k1*L1);
Z(6,2)      = -beta1*i*k1^3*exp(i*k1*L1);
Z(6,3)      = -beta1*k1^3*exp(-k1*L1);
Z(6,4)      = beta1*i*k1^3*exp(-i*k1*L1);
Z(6,5)      = -beta2*k2^3;
Z(6,6)      = beta2*i*k2^3;
Z(6,7)      = beta2*k2^3;
Z(6,8)      = -beta2*i*k2^3;

Z(7,1:4)    = [0 0 0 0];
Z(7,5)      = beta2 * k2^2 * exp(k2 * L2);
Z(7,6)      = -beta2 * k2^2 * exp(i * k2 * L2);
Z(7,7)      = beta2 * k2^2 * exp(-k2 * L2);
Z(7,8)      = -beta2 * k2^2 * exp(-k2 * i * L2);

Z(8,1:4)    = [0 0 0 0];
Z(8,5)      = beta2 * k2^3 * exp(k2 * L2);
Z(8,6)      = -beta2 * i * k2^3 * exp(i * k2 * L2);
Z(8,7)      = -beta2 * k2^3 * exp(-k2 * L2);
Z(8,8)      = beta2 * i * k2^3 * exp(-i * k2 * L2);           

E1          = exp( -(1+i)*L1 / delta );
E2          = exp( -(1+i)*L2 / delta );

u(1,1)      = -(s1 + t1);
u(2,1)      = -c*t1 + alph*c^2*t1;
u(3,1)      = -(s1 + t1*E1) + (s2 + t2);
u(4,1)      = -(c*t1*E1) + (c*t2);
u(5,1)      = -beta1*c^2*t1*E1 + beta2*c^2*t2;
u(6,1)      = -beta1*c^3*t1*E1 + beta2*c^3*t2;
u(7,1)      = -beta2*c^2*t2*E2 + phi;
u(8,1)      = -beta2*c^3*t2*E2 + zta;

r = Z\u;

%r = inv(Z)*u;


