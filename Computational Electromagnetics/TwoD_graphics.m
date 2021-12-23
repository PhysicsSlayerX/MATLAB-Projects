xa = linspace(-2,2,50);
ya = linspace(-1,1,25);
[Y,X] = meshgrid(ya,xa);
D = X.^2 + Y.^2;
h = imagesc(xa,ya,D');
h2 = get(h,'Parent');
set(h2, 'YDir', 'normal');
axis equal tight;