xa = linspace(-2,2,50);
ya = linspace(-1,1,25);
[Y,X] = meshgrid(ya,xa);

D = X.*sin(X.^2) - Y.*cos(Y.^3);
surf(xa,ya,D');

axis equal tight;
shading interp;

camlight; lighting phong;
view(25,45);
%colormap('jet(1024)');

CMAP = zeros(256,3);
c1 = [0 0 1]; %blue
c2 = [1 1 1]; %white
c3 = [1 0 0]; %red
for nc = 1 : 128
    f = (nc - 1)/128;
    c = (1 - sqrt(f))*c1 + sqrt(f)*c2;
    CMAP(nc,:) = c;
    c = (1 - f^2)*c2 + f^2*c3;
    CMAP(128+nc,:) = c;
end
colormap(CMAP);