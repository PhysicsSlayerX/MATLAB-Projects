%Triple Scalar Product

degrees = pi/180;
A = [ 45*cos(10*degrees) -45*sin(10*degrees)];
B = [ -67*cos(20*degrees) 67*sin(20*degrees)];
C = [ -23*cos(69*degrees) -23*sin(69*degrees)];
D = [ -98*cos(32*degrees) -98*sin(32*degrees)];

%D = cross(A, cross(B,C));

E = A - B + C - D;
Emag = sqrt(E(1)^2 + E(2)^2);
theta = atan(E(2)/E(1))./degrees;


A = [1  1  1];
B = [1 -1 -1];
C = [-1 1 -1];
D = [1  1 -2];
E = [2  2  2];

Find = cross(A, A)