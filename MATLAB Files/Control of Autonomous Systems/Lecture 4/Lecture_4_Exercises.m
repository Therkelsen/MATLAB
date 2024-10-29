clc; clear; close all;
format compact

disp('State Space Model')
A = [-0.4, 0, -0.01;
     1, 0, 0;
     -1.4, 9.8, -0.02]F
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

% Desired poles for the closed-loop system
disp(' ')
disp(' ')
disp('State Feedback Control')
lambda_cl = [-6, -3 + 2i, -3 - 2i]

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

u_star = 0;

% reference
x_star = [0;   % q(t): pitch angular velocity (rad/s)
          0;   % theta(t): pitch angle (rad)
          10]; % u(t): longitudinal velocity (m/s)

% initial cond (same units)
init_cond = [0.1, pi/8, 0];

disp(' ')
disp(' ')
disp('LQR Control')
% Define Q and R matrices
% A common choice is to penalize all states equally
Q = C' * C
% Scalar R, penalizes the control effort
R = 1

% Compute the LQR gain matrix K
K_lqr = lqr(A, B, Q, R)

disp(' ')
disp(' ')
disp('Discretized LQR Control')

T = 0.01;
sample_period = 10*T;

sys_d = c2d(sys, sample_period, 'zoh')

% Extract the discrete-time matrices
[A_D, B_D, C_D, D_D] = ssdata(sys_d);

% Compute the discretized LQR gain matrix K
K_lqr_d = dlqr(A_D, B_D, Q, R)

disp(' ')
disp(' ')
disp('LQI Control with Disturbance Rejection - Revised')

















