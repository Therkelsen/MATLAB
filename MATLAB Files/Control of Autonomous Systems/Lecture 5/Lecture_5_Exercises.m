clc; clear; close all;
format compact
% sampling time
Ts = 0.01;

% State-space representation
A = [ -0.4 ,   0 , -0.01  ;
        1  ,   0 ,     0  ;
      -1.4 , 9.8 , -0.02  ]; 

B = [ 6.3 ; 0 ; 9.8 ];

p = eig(A) % eigenvalues of system (in open-loop)

C = [ 0 , 0 , 1 ];

D = 0;

% initial cond (same units)
init_cond = [0.1, pi/8, 0];

sysCT = ss(A,B,C,D); % creates a system object in CT

% Linear Quadratic Regulator
Q = 10*eye(3);
R = eye(1);
K = lqr(A,B,Q,R);

% 1.3
% Linear Quadratic Regulator with integration
Q_E = 10*diag([1;1;1;1]);
R = eye(1);
K_E = lqi(sysCT,Q_E,R)

% 2.1
T = 0.01;
sample_period = 10*T;

sys_d = c2d(sysCT, sample_period, 'zoh')

% Extract the discrete-time matrices
[A_D, B_D, C_D, D_D] = ssdata(sys_d);

Q_E_D = 10*diag([1;1;1;1]);
R_D = eye(1);
K_E_D = lqi(sys_d, Q_E_D, R_D);