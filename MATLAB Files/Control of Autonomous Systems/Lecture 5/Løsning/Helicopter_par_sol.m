% Helicopter parameters
% Jerome Jouffroy, September 2024

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

sysCT = ss(A,B,C,D); % creates a system object in CT



% Linear Quadratic Regulator
Q = 10*eye(3);
R = eye(1);
K = lqr(A,B,Q,R);

% Linear Quadratic Regulator with integration
Q_E = 10*diag([1;1;1;1]);
R = eye(1);
K_E = lqi(sysCT,Q_E,R);


