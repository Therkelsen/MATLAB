clc; clear; close all;
format compact

utils = Utils;

disp("Problem 2.1")

% The file ‘dataset_problem_2_1.mat’ contains the data matrix
% for a random sample of 5-dimensional vector observations
% from a population possibly being multivariate normal. 

% Perform a model check for multivariate normality:
%   Calculate the descriptive statistics x^bar, S, R
%   Make plots of the multivariate scatter matrix plus marginal
%      histograms and boxplots
%   Calculate the squared sample Mahalanobis distances,
%      d^2_i = (x_i - x^bar)^T * S^-1 * (x_i - x^bar), and make
%      a QQ-plot versus the relevant X^2 distribution for check
%      of approximate linearity
%   Conclude

data = load("Lecture 2/dataset_problem_2_1.mat").X;
fprintf('Data Size: [%d %d]\n', size(data, 1), size(data, 2));

[mean, S, range] = utils.calculate_descriptive_statistics(data, true, false);

mahalanobis_distances = utils.calculate_mahalanobis_distances(data, mean, S, true, false);

% The dataset does not appear linear like the chi squared distribution. It
% must be a different kind of distribution, like an exponential one.
%%
clc; clear; close all;
format compact

utils = Utils;

disp("Problem 2.2")

% Perform the eigenvalue decomposition of Σ, i.e. Σ P Λ P^T

mu = [1;
      2;
      0]

Sigma = [4, 1, 2;
         1, 9, -3;
         2, -3, 5]

[P, lambda] = eig(Sigma)

% Find Σ^(1⁄2) = P Λ^(1⁄2) P^T
% Where P = matrix of eigen vectors, Λ is diagonal matrix of eigenvalues

Sigma_sqrt = P * sqrt(lambda) * P'

n = 10000;
p = 3;
Y = randn(n, p);

% Generate Xi = μ + Σ^(1⁄2)*Yi, be sure to keep track of dimensions
X = mu' + Y*Sigma_sqrt;

% Calculate the descriptive statistics
% Make plots of the multivariate scatter matrix plus marginal histograms and boxplots
[mu_hat, Sigma_hat, range] = utils.calculate_descriptive_statistics(X, true, false);

% Calculate the squared sample Mahalanobis distances,
% and make a QQ-plot versus the relevant distribution
% for check of approximate linearity

mahalanobis_distances = utils.calculate_mahalanobis_distances(X, mu_hat, Sigma_hat, true, false);

% Conclusion: sometimes maybe good sometimes maybe shit. this time shit.