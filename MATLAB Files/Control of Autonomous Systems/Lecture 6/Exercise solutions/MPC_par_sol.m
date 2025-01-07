% DT double-integrator MPC example (for Simulink implementation)
% Jerome Jouffroy, October 2023

A = [ 1 1 ; 0 1 ];
B = [ 0 ; 1 ];
C = [ 1 0 ];
x0 = [ 10 ; 0 ];

Q = C'*C;
R = 1/10;

Qb = blkdiag(Q,Q,Q);
Rb = blkdiag(R,R,R);

Umin = -1*ones(3,1);
Umax = 1*ones(3,1);

Ab = [ eye(2) ; A ; A^2 ];
Bb = [ zeros(2,3) ; B zeros(2,2) ; A*B B zeros(2,1) ];

H = 2*(Bb'*Qb*Bb + Rb);
preF = Ab'*Qb*Bb;

N = 41; % number of iterations

