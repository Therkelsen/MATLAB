clc, clear

% Constants a = [q1, q2, q3, J]

q1 = 1;
q2 = 10;
q3 = 0.1;
J = 40;

a = [q1, q2, q3, J];

x0 = [2; 200];

% EQ point:

x2_star = 300;
x1_star = q3/q2 * x2_star;
u_star = q1 * x1_star * x2_star;

x_star = [x1_star; x2_star];

% Linearization matrices

A = [-q1 * x2_star, -q1 * x1_star; 
     q2/J, -q3/J];

B = [1; 0]
C = [1, 0];
D = 0;


K = place(A, B, [-10, -11]);
