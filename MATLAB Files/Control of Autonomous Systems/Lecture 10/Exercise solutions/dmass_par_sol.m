% Double-mass parameters and motion planning
% Jerome Jouffroy, November 2024

m1 = 10;
m2 = 20;
k = 5;

% state-space representation
A = [ 0 0 1 0 ; 0 0 0 1 ; -k/m1 k/m1 0 0 ; k/m2 -k/m2 0 0 ];
B = [ 0 ; 0 ; 1/m1 ; 0 ];

% controllability study
Wc = [ B , A*B , A^2*B , A^3*B ];
iscontrollable = det(Wc);

Bt = [ 0 ; 0 ; 0 ; 1 ];

% finding the T transform
T1 = Bt'*inv(Wc)
T2 = T1*A;
T3 = T2*A;
T4 = T3*A;

T = [ T1 ; T2 ; T3 ; T4 ]

At = T*A*inv(T)

% motion planning / feedforward / open-loop control part
% we want to steer both masses to a position of 10, starting from 0, in 10
% seconds
x0 = [ 0 ; 0 ; 0 ; 0 ];
xT = [ 40 ; 40 ; 0 ; 0 ];
H = 10;

z0 = T*x0;
zT = T*xT;

% calculation of polynomial coefficients
alpha_vect = poly7traj(z0,zT,H);

% extraction of system coefficients (in the z-space)
a_vect = -At(4,1:4)';
