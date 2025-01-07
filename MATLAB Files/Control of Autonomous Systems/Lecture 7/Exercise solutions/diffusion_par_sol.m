% parameter for the diffusion observer
% Jerome Jouffroy, November 2024


% state-space representation modelling
dim = 10;
A = -2*eye(dim);
A(2:dim,1:dim-1) = A(2:dim,1:dim-1) + eye(dim-1);
A(1:dim-1,2:dim) = A(1:dim-1,2:dim) + eye(dim-1);
A(1,1) = -1; % Neumann conditions
A(dim,dim) = -1;

C = zeros(1,dim);
C(1) = 1;


% look at the stability of matrix A
lambda = eig(A)';



% calculation of different sets of observer gains
p_obs_brute = -1*ones(dim,1);
L_brute= acker(A',C',p_obs_brute)';
lambda_brute = eig(A-L_brute*C)';

L_naive = zeros(dim,1);
L_naive(1) = 10;
lambda_naive = eig(A-L_naive*C)';

lambda = eig(A);
lambda_des = lambda;
lambda_des(dim) = -2;
L_des= acker(A',C',lambda_des)';
lambda_des_obs = eig(A-L_des*C)';

