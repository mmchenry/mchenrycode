function visModels(action,freq)
% Runs the numerical model of the neuromast


if nargin < 1
    action = 'animate comparison';
end

cupRoot = ['/Network/Servers/atlantis.local/Users/mjm/Documents' ...
            '/Projects/nmast_analytical_model'];

switch action
    case 'series'
        visModels('animate basic',1);
        visModels('animate basic',44);
        visModels('animate basic',100);
        
    case 'animate basic'
        cd([cupRoot filesep 'movies'])

        c           = c_default_theo;  
        
        if nargin <2
            c.freqs     = 1;
        else
            c.freqs     = freq;
        end
        
        c       = numerical_twopart(c,'torsion spring');
        time    = linspace( 0 , 1/c.freqs , 60 );
        
        fName   = ['frq_' num2str(c.freqs) '_Hz'];
        
        aviobj = avifile([cupRoot filesep 'movies' filesep fName '.avi'],...
            'Compression','None','Quality',50,'fps',5);
        
        y           = c.heights;
        x           = c.M;
        omega       = 2.*pi.*c.freqs;
        delta       = (2*c.mu/(omega*c.rho))^0.5;

        xMax        = 1.1*max(abs(x));
        yMax        = 1.0*max(y);
        numCol      = 7;
        numRow      = 10;
        vectX       = -xMax:xMax/numCol:xMax;
        vectY       = yMax/numRow:yMax/numRow:yMax;
        
        flwAmp      = xMax.*i.*omega;   
        bLayer      = flwAmp .* (1-exp(-y.*(1+i)/delta)) ;
        
        bLayer_norm = bLayer / abs(bLayer(end));
        
        F           = figure('DoubleBuffer','on');

        for j = 1:length(time)
            xCup        = real(x .* exp(i*omega*time(j)));
            xFlow       = real(bLayer_norm .* exp(i*omega*time(j)));
            
            uWater      = interp1(y,xFlow,vectY);
            vWater      = eps.*ones(size(vectX));

            [xMesh,yMesh] = meshgrid(vectX,vectY);
            [vMesh,uMesh] = meshgrid(vWater,uWater);
            
            sclFactor   = .1 * 2*xMax;
            
            quiverc(xMesh,yMesh,uMesh,vMesh,'AbsScale',...
                sclFactor,'w-','HeadSize',.4);hold on;

            set(gca,'Color','k')
            set(gcf,'Color','w')
            
            clr = [0.3 1 0.3];
            h2  = plot(xCup,y);hold off
            set(h2,'LineWidth',7,'Color',clr)
            axis([-xMax xMax 0 1.1*yMax]);
            title(['freq = ' num2str(c.freqs) ' Hz']);
            frame       = getframe(F);
            aviobj      = addframe(aviobj,frame);

            pause(.01);
            hold off
        end
        aviobj = close(aviobj);
        close
        disp('done saving')
        
      
    case 'animate comparison'
        load('theoData')
        theoXs      = [ theo.pt1.xAvg  theo.pt2.xAvg     theo.pt3.xAvg     ];
        theoY       = [ theo.pt1.y     theo.pt2.y        theo.pt3.y        ];
        time        = theo.time.avg;


        cd([cupRoot filesep 'movies'])

        c           = c_default_theo;
        
        c.kinoHeight= 16e-6;
        c.dispAmp       = .26 * 10^-6; %m
        %c.dispAmp       = .23 * 10^-6; %m
        c.freqs         = 44;
        c.numHairs      = 11;

        c.cupHeight     = 45e-6;
        

        c       = numerical_twopart(c,'torsion spring');
        %time    = linspace( 0 , 1/c.freqs , 60 );

        fName   = ['comparison_' num2str(c.freqs) '_Hz'];

        aviobj = avifile([cupRoot filesep 'movies' filesep fName '.avi'],...
            'Compression','None','Quality',100,'fps',5);

        y           = c.heights;
        x           = c.M;
        xMax        = 1.1*max([abs(x);max(theoXs(:))]);
        yMax        = 1.0*max(y);
        omega       = 2.*pi.*c.freqs;

        F           = figure('DoubleBuffer','on');

        for j = 1:length(time)
            xCup        = real(x .* exp(i*omega*time(j)));
            
            clr = 0.5.*[1 1 1];
            h2  = plot(xCup,y);hold on
            h1 = plot(theoXs(j,:),theoY,'ro-');
            hold off
            set(h2,'LineWidth',7,'Color',clr)
            axis([-xMax xMax 0 1.1*yMax]);
            title(['freq = ' num2str(c.freqs) ' Hz']);
            frame       = getframe(F);
            aviobj      = addframe(aviobj,frame);

            pause(.01);
            hold off
        end
        aviobj = close(aviobj);
        close
        disp('done saving')
end



function c = c_default_theo
%Parameters for all anlayses
 c.numHeights    = 50;  
 c.bunHeight     = 5.3e-6; %From Dinklo, 2005
 c.dispAmp       = 10 * 10^-6; %m
 c.E_matrix      = 31; %31 Pa
 c.EI_kino       = 2e-21; % 2e-21 N m^2
 c.bundleStiff   = 2.925e-14; %Nm/rad (van Netten & Kroese, 1987) 
 c.linStiff      = 0.13 * 10^-3; %N/m (van Netten & Kroese, 1987) 
 c.rho           = 998; %998 kg m^-3
 c.mu            = 1.002e-3; %1.002e-3 Pa s

%Data from morphometric measurements (based on stiffness paper)
c.baseDiameter 	= 8.88e-6;
c.midDiameter 	= 7.2e-6 ;
%c.kinoHeight 	= 29.7e-6;
c.kinoHeight 	= 16e-6;
c.cupHeight     = 45e-6;
c.numHairs      = 11;


function c = runAnalytical(c,mType)

%Define general parameters
E_matrix    = c.E_matrix;
EI_kino     = c.EI_kino;
bundleStiff = c.bundleStiff;
rho         = c.rho;
mu          = c.mu; 
D0          = c.dispAmp;
L           = c.cupHeight;
a           = mean([c.baseDiameter c.midDiameter])/2;
freqs       = c.freqs;

I           = (pi / 64) * (2 * a)^4;
EI          = (I * E_matrix) + (EI_kino * c.numHairs);
baseStiff   = bundleStiff .* c.numHairs;
linStiff    = c.linStiff .* c.numHairs;

num1        = round(c.numHeights .* c.kinoHeight / c.cupHeight);
x1          = [linspace(0,c.kinoHeight, num1)]';
x1          = x1(2:end);
x2          = [linspace(c.kinoHeight,c.cupHeight,c.numHeights-num1+1)]';
x2          = x2 - x2(1);
c.heights   = [x1(1:end-1);x2+x1(end)];

for ii = 1:length(freqs)
    omega       = 2*pi*freqs(ii);
    if strcmp(mType,'low')
        c.M(:,ii)   = simpleBeam_low_2pin(c.heights,omega,mu,rho,L,D0...
                                        ,a,EI,baseStiff);
    elseif strcmp(mType,'high') 
        c.M(:,ii)   = simpleBeam_high(c.heights,omega,mu,rho,L,D0...
                                        ,a,EI,linStiff);
    elseif strcmp(mType,'bLayer') 
        c.M(:,ii)   = bLayer(c.heights,omega,mu,rho,D0);
    else
        error('incorrect model type');
    end
end

function v = simpleBeam_low(x,omega,mu,rho,L,D0,a,EI,kn)
%one-part model for low freq - all fluid forces

delta   = (2 * mu / rho / omega)^0.5;
s       = a ./ 2^.5 / delta;
g       = 0.577 + log(s);
G       = -g / (g^2 + (pi/4)^2);
d       = D0 * (1 + i) * x ./ delta;
b       = omega * pi * (2*rho*a^2*omega - 4*mu*G*i - pi*mu*G/g);
%b       = omega * pi * (- 4*mu*G*i);
v       = -(1/120 + i/120) * b*D0.*x .* ...
            (40*EI*L^3 + kn.*x .* (20*L^3 - 10*L^2.*x + x.^3)) ./ ...
            (EI * kn *delta);

function v = simpleBeam_low_2pin(x,omega,mu,rho,L,D0,a,EI,kn)
%one-part model for low freq - all fluid forces - anchored base

delta   = (2 * mu / rho / omega)^0.5;
s       = a ./ 2^.5 / delta;
g       = 0.577 + log(s);
G       = -g / (g^2 + (pi/4)^2);
d       = D0 * (1 + i) * x ./ delta;
b       = omega * pi * (2*rho*a^2*omega - 4*mu*G*i - pi*mu*G/g);
v       = -(1/120 + i/120) * b*D0.*x.^2 .* ...
            (20*L^3-10*L^2.*x+x.^3) ./ (EI*delta);        
        
% function v = simpleBeam_high(x,omega,mu,rho,L,D0,a,EI,kn)
% %one-part model for low freq - all fluid forces
% 
% delta   = (2 * mu / rho / omega)^0.5;
% s       = a ./ 2^.5 / delta;
% g       = 0.577 + log(s);
% G       = -g / (g^2 + (pi/4)^2);
% d       = D0 * (1 + i) * x ./ delta;
% b       = omega * pi * (2*rho*a^2*omega - 4*mu*G*i - pi*mu*G/g);
% %b       = omega * pi * (2*rho*a^2*omega - pi*mu*G/g);
% v       = - b*D0.*x .* (12*EI*L^2 + kn.*x .* (6*L^2 - 4*L.*x + x.^2)) ./ ...
%             (24.*EI.*kn);

function v = simpleBeam_high(x,omega,mu,rho,L,D0,a,EI,k)
%one-part model for low freq - all fluid forces

delta   = (2 * mu / rho / omega)^0.5;
s       = a ./ 2^.5 / delta;
g       = 0.577 + log(s);
G       = -g / (g^2 + (pi/4)^2);
d       = D0 * (1 + i) * x ./ delta;
%b       = omega * pi * (2*rho*a^2*omega - 4*mu*G*i - pi*mu*G/g);
b       = omega * pi * (2*rho*a^2*omega - pi*mu*G/g);
v       = - b*D0.*L .*(6.*EI + k.*x.^2 .* (-3.*L + x)) ./ ...
            (6*EI*k);

function U = bLayer(x,omega,mu,rho,D0)

delta   = (2*mu./(rho.*omega)).^0.5;
disp    = D0 .* (1 - exp(-x.*(1+i) ./ delta));
U       = disp.*i.*omega;









