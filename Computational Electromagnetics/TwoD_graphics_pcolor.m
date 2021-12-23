xa = linspace(-2,2,50);
ya = linspace(-1,1,25);
[Y,X] = meshgrid(ya,xa);
D = X.^2 + Y.^2;
pcolor(xa,ya,D');
shading interp;
axis equal tight;