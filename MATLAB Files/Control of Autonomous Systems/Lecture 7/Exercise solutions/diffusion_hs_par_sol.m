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

B = zeros(dim,1);

% adding a head source
ihs = 7; % index heat source
B(ihs) = 1;
A(ihs,ihs) = A(ihs,ihs) -1;

% look at the stability of matrix A
eig(A);

% calculation of different sets of observer gains
p_obs_brute = -1*ones(dim,1);
L_brute= acker(A',C',p_obs_brute)';

L_naive = zeros(dim,1);
L_naive(1) = 10;

p_ol = eig(A);
p_obs_des = p_ol;
p_obs_des(dim) = -2;
L_des= acker(A',C',p_obs_des)';
