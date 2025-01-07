% Pendulum parameters
% Jerome Jouffroy, September 2023

m = 50;
l = 1;
g = 9.8;
d = 20;

% parameter vector
a = [ m ; g ; l ; d ];

% tuning for stabilizing part
A = [ 0 , 1 ; 0 , 0 ]; % ie we have double-integrator dynamics
B = [ 0 ; 1 ];
p_des = [ -1 ; -1 ];
K = acker(A,B,p_des);

