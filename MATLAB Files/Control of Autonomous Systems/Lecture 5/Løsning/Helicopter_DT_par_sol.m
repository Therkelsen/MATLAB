% Helicopter parameters
% Jerome Jouffroy, October 2023

A = [ -0.4 ,   0 , -0.01  ;
        1  ,   0 ,     0  ;
      -1.4 , 9.8 , -0.02  ]; 

B = [ 6.3 ; 0 ; 9.8 ];

C = [ 0 , 0 , 1 ];

D = 0;

sys_CT = ss(A,B,C,D);

% sampling period
Ts = 0.01;

% discretization
sys_DT = c2d(sys_CT,Ts);

% Linear Quadratic Regulator with integration (discrete-time)
Q = 10*eye(4);
R = eye(1);
K = lqi(sys_DT,Q,R);

