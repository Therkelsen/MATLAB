% Problem 1.1

clc; clear; close all;

disp('Problem 1.1')

% Given covariance matrix
cov = [4, 1, 2;
       1, 9, -3;
       2, -3, 25]

% Using eig
disp('Problem 1.1.1')
[eig_vec_, eig_val_] = eig(cov)

% Confirm that cov*eig_vec = eig_vec*eig_val
cov*eig_vec_ - eig_vec_*eig_val_

% Doing it manually
syms l
I = eye(size(cov))
char_eq = det(cov - l * I)
% Eigen values
eig_val = diag(double(solve(char_eq, l)))

% Eigen vectors
cov_lambda1 = cov - eig_val(1,1) * I;
eig_vec_1 = null(cov_lambda1);
eig_vec_1 = eig_vec_1 / norm(eig_vec_1);

cov_lambda2 = cov - eig_val(2,2) * I;
eig_vec_2 = null(cov_lambda2);
eig_vec_2 = eig_vec_2 / norm(eig_vec_2);

cov_lambda3 = cov - eig_val(3,3) * I;
eig_vec_3 = null(cov_lambda3);
eig_vec_3 = eig_vec_3 / norm(eig_vec_3);

eig_vec = [eig_vec_1, eig_vec_2, eig_vec_3]

disp('Problem 1.1.2')
% Verify that |Sigma| = Pi^3_i=1 (lambda_i)
eig_val_prod = prod(eig_val);
eig_val_prod_ = prod(eig_val_);
det(cov);

det(cov) - eig_val_prod
det(cov) - eig_val_prod_

% False?

disp('Problem 1.1.3')
% Verify that trace(cov) = Pi^3_i=1 (lambda_i)
trace(cov) - eig_val_prod
trace(cov) - eig_val_prod_

% Use EVD and square root to find cov^(1/2) and cov^(-1/2)
D_sqrt = diag(sqrt(eig_val));
D_inv_sqrt = diag(1 ./ sqrt(eig_val));

% Compute cov^(1/2)
cov_sqrt = eig_vec .* D_sqrt .* eig_vec'

% Compute cov^(-1/2)
cov_inv_sqrt = eig_vec .* D_inv_sqrt .* eig_vec'
%% 

% Problem 1.2

clc; clear; close all;

utils = Utils;

disp('Problem 1.2')

% Loading data which contains track records for women in 54 countries for
% the running disciplines 100m [s], 200m [s], 400m [s], 800m [min], 1500m
% [min], 3000m [min], marathon [min].
 
data = load("Lecture 1\dataset_problem_1_2.dat");

% I'm gonna convert the second columns into minutes, for the graphs.
data(:, 1:3) = data(:, 1:3) / 60;

% Original headers and data
headers = {'100m [min]', '200m [min]', '400m [min]', '800m [min]', '1500m [min]', '3000m [min]', 'marathon [min]'};
data = array2table(data, 'VariableNames', headers)

disp('Problem 1.2.1')

descriptive_statistics = utils.calculate_descriptive_statistics(table2array(data), headers);

% Looking for ellipses determining proportionality between the different
% data. For example box 2 from top left is x1 as a function of x2. Other
% than a few outliers, this has very good proportionality, meaning that
% if you are good at 200m running, you are likely good at 100m running too.
