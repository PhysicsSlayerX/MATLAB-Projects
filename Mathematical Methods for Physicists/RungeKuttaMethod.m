% ETD1RK for a simple bubble shrinking Allen-Cahn dynamics
%    phi_t = eps*lap(phi) - 1/eps*W'(phi)
%Parameter Values
epsilon = 0.05;
epsInv  = 1/epsilon;
kappa   = 300;

Nmod  = 100;
dt    = 5e-3;
T     = 100;
%Initial set-up
m      = 2^8;
n      = 2^8;
Lx     = 1.5;
Ly     = 1.5;
x      = linspace(-Lx, Lx, m+1); x = x(1:end-1); dx = 2*Lx/m;
y      = linspace(-Ly, Ly, n+1); y = y(1:end-1); dy = 2*Ly/n;
[x,y]  = meshgrid(x,y);
%Initial phi
% Circular Initial
r0 = 1; a = 1; b = 1;
r  = sqrt(x.^2/a^2+y.^2/b^2);
phi0 = 0.5 + 0.5*tanh(3*(r0-r)/epsilon);
%Fourier Spectral Setting
k1 = [0:m/2, -(m/2-1):-1];
k2 = [0:n/2, -(n/2-1):-1];
[k1, k2] = meshgrid(k1, k2);
k1Sqr  = k1.^2;
k2Sqr  = k2.^2;
hx = (pi/Lx)^2;
hy = (pi/Ly)^2;
%Other Pre-calculated Variables
lambdaij = -k1Sqr.*hx - k2Sqr.*hy;
lij = epsilon*lambdaij - kappa;
lijInv = 1./lij;
expTerm = exp(lij*dt);
%For Loop Iteration
for i = 1 : T/dt

     G  = 18*(phi0-phi0.^2).^2;
     Gp = 36*(phi0-phi0.^2).*(1-2*phi0);
     phiHat0 = fft2(phi0);

     Nonlin    = kappa*phi0 - epsInv*Gp;
     NonlinHat = fft2(Nonlin);

     phiHat = expTerm.*phiHat0 + lijInv.*(expTerm-1).*NonlinHat;

     phi  = ifft2(phiHat);
     phi0 = phi;

     if (mod(i, Nmod)==0)
         figure(1)
         drawnow;
         pcolor(x,y,phi); shading('interp');
         axis image; axis off;
         colorbar('SouthOutside');
         colormap(jet);

         time = i*dt;
         fprintf('Time = %1.4f\v eps = %1.4f\v kappa = %1.4f\n',...
                  time,epsilon,kappa);
     end

end

