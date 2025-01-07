% Parameters for the Mass-Spring-Damper
% as well as parameters of the State-Variable Filters
% Jerome Jouffroy, November 2023

m = 10;
k = 5;
d = 0.1;

% definition of associated parameters
a0 = k/m;
a1 = d/m;
b = 1/m;

% parameter vector (ie the one to be estimated)
theta = [ a0 ; a1 ; b ]

% definition of state-variable filters
A = [ 0 1 0 ; 0 0 1 ; 0 0 0 ];
B = [ 0 ; 0 ; 1 ];
poles = -10*ones(3,1);
K = acker(A,B,poles);
A_F = A-B*K;

C_F0 = [ 1 0 0 ];
C_F1 = [ 0 1 0 ];
C_F2 = [ 0 0 1 ];

B_F = [ 0 ; 0 ; -A_F(3,1)]