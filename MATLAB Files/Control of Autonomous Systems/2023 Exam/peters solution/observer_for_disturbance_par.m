clc, clear

% Constants a = [q1, q2, q3, J]

q1 = 1;
q2 = 10;
q3 = 0.1;
J = 40;

a = [q1, q2, q3, J];

x0 = [2; 200; 0];

% First write statespace since q1 is 0
% Define new state: x = [m, w, d]

A = [0, 0; q2/J, -q3/J];
B = [1; 0];
C = eye(2,2);
D = zeros(1, 2);

Ae = [A, B; zeros(1,3)]
Be = [B; 0]
Ce = [C, zeros(2,1)]


lambda_des = [-10, -12, -8]
L_des= place(Ae',Ce',lambda_des)'

lambda_des_obs = eig(Ae-L_des*Ce)';