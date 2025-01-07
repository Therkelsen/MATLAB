% parameters for the satellite antenna problem
% Jerome Jouffroy, November 2024

% system parameters
a_m = 1.71;
Kl = 100;
Kg = 0.1;

% state-space representation
A = [ 0 1 ; 0 -a_m ];
B = [ 0 ; Kl*Kg ];
C = [ 1 , 0 ];
D = 0;

% stability of original system
p_sys = eig(A)

% controller tuning
p_cont = p_sys(2)*ones(2,1);
K = acker(A,B,p_cont);

% feedforward gain computation
N = inv([ A , B ; C , 0 ]);
N_x = N(1:2,3);
N_u = N(3,3);
N_bar = N_u + K*N_x;

% observer tuning
p_obs = 3*p_cont;
L = acker(A',C',p_obs)';