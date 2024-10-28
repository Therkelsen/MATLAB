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



% Linear Quadratic Regulator (with different tunings)
Q = eye(3);
R = 100*eye(1);
%R = 1000*eye(1); % different tuning
%R = 1e5*eye(1);
K = lqr(A,B,Q,R);

% feedforward gain computation
N = inv([ A , B ; C , 0 ]);
N_x = N(1:3,4);
N_u = N(4,4);

N_bar = N_u + K*N_x;

% discretization
sysDT = c2d(sysCT,10*Ts,'tustin');
[Ad,Bd,Cd,Dd] = ssdata(sysDT);

% stability of DT system in OL?
lambda_DT = eig(Ad) % not stable because 2 of the eigenvalues have their norm bigger than 1.

% Linear Quadratic Regulator in discrete-time
Kd = dlqr(Ad,Bd,Q,R);

% feedforward gain computation in discrete-time
Nd = inv([ Ad-eye(3) , Bd ; Cd , Dd ]);
Nd_x = Nd(1:3,4);
Nd_u = Nd(4,4);

Nd_bar = Nd_u + Kd*Nd_x;
