% Pendulum parameters
% Jerome Jouffroy, September 2024

m = 50;
l = 1;
g = 9.8;
d = 0.1;

% parameter vector
a = [ m ; g ; l ; d ];


% equilibrium points

xstar = [ pi/4 ; 0 ];

ustar = m*g*l*sin(xstar(1));

% initial condition
x0 = xstar;

% linear approximation of pendulum around xstar (and associated ustar)
A = [ 0 1 ; -g/l*cos(pi/4) -d/(m*l^2) ];
B = [ 0 ; 1/(m*l^2) ];