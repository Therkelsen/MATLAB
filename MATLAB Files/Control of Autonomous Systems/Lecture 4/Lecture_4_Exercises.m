clc; clear; close all;
format compact

% Exercise 1, 2
A = [-0.4, 0, -0.01;
     1, 0, 0;
     -1.4, 9.8, -0.02]

B = [6.3;
     0;
     9.8]

% Since C (output) is not explicitly stated, we can assume
% that it is an identity matrix same dimensions as A
C = eye(size(A))

% Feedthrough is not defined, so this is a zero vector
% with the same size as B
D = zeros(size(B))

sys = ss(A, B, C, D);

e = eig(A)
poles = real(e);
% Check the real parts of the eigenvalues
disp(" ")
disp('Real parts of the eigenvalues: ')
disp(poles)
disp(" ")
if all(poles < 0)
    disp('The system is stable.');
elseif any(poles > 0)
    disp('The system is unstable.');
else
    disp('The system is marginally stable.');
end

% Exercise 3 solved in simulink

% Exercise 4

% Desired poles for the closed-loop system
lambda_cl = [-2, -1 + 1i, -1 - 1i]

K = place(A, B, lambda_cl)

A_cl = A - B * K

e_cl = eig(A_cl)
poles_cl = real(e_cl);
% Check the real parts of the eigenvalues
disp(" ")
disp('Real parts of the eigenvalues: ')
disp(poles_cl)
disp(" ")
if all(poles_cl < 0)
    disp('The system is stable.');
elseif any(poles_cl > 0)
    disp('The system is unstable.');
else
    disp('The system is marginally stable.');
end