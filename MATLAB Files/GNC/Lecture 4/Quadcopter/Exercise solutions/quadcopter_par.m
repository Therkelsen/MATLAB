% Quadcopter parameters, taken from Lecture 2, question 1
% Jerome Jouffroy, February 2024

m = 4; % mass in kg
g = 9.81; % gravity constant

% inertia matrix
Ic = diag([ 2 ; 2 ; 4 ]);

% rigid-body mass matrix
M_RB = [ m*eye(3) , zeros(3) ; zeros(3) , Ic ];