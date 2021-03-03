syms a
syms b
syms g
a = deg2rad(30);
b = deg2rad(30);
g = deg2rad(30);

syms Rz

syms Ry

syms Rx

Rz = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1]

Ry = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)]

Rx = [1 0 0; 0 cos(g) -sin(g); 0 sin(g) cos(g)]
 
RzRyRx = Rz*Ry*Rx
