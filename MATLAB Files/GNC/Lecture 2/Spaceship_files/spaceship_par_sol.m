% Spaceship parameters, taken from Crew Dragon
% Jerome Jouffroy, February 2022

m = 12055; % mass in kg
H = 8.1; % height in m (with trunk)
R = 2; % radius base in m

% compute inertia matrix, approximating spacecraft as cylinder
Ic = diag([ 1/12*m*H^2+1/4*m*R^2 ; 1/12*m*H^2+1/4*m*R^2 ; 1/2*m*R^2 ]);

% rigid-body mass matrix
M_RB = [ m*eye(3) , zeros(3) ; zeros(3) , Ic ];