% parameter file for the ADRC pendulum simulation
% Jerome Jouffroy, October 2023

m = 50;
l = 1;
g = 9.8;
d = 0.1;

% parameter vector
a = [ m ; g ; l ; d ];

% ADRC parameters
b0 = 1/m; % because l=1
T_settle = 4;
pole_CL = 6/T_settle;
KP = pole_CL^2;
KD = 2*pole_CL;

pole_ESO = -15;
poles_ESO = pole_ESO*ones(3,1);

% observer matrices
A = [ 0 , 1 , 0 ; 0 , 0 , 1 ; 0 , 0 , 0 ];
B = [ 0; b0 ; 0 ];
C = [ 1 , 0 , 0 ];
L = acker(A',C',poles_ESO)';