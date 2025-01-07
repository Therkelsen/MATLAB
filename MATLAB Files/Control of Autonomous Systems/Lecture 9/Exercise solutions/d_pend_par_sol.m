% Pendulum parameters
% Jerome Jouffroy, November 2024

m = 50;
l = 1;
g = 9.8;
d = 0.1;

% parameter vector (for the plant model)
a = [ m ; g ; l ; d ];

% parameter definitions for state estimation
a0 = g / l;
a1 = d / (m*l^2);
b = 1 / (m*l^2);