e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];

e1p = [0.33;0;0];
e2p = [0;0.33;0];
e3p = [0;0;0.33];

A1 = 10;
A2 = 5;
A3 = 3;

A = A1*e1 + A2*e2 +A3*e3;
A = [A1;A2;A3];
Ap = [dot(e1,e1p) dot(e2,e1p) dot(e3,e1p);
      dot(e1,e2p) dot(e2,e2p) dot(e3,e2p);
      dot(e1,e3p) dot(e2,e3p) dot(e3,e3p)]*A; 