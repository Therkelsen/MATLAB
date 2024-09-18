clc; clear; close all;

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

data = load("Lecture 2\dataset_problem_2_1.mat").X;

% Define column names based on the number of columns in your data
headers = arrayfun(@(x) sprintf('Dim %d', x), 1:size(data, 2), 'UniformOutput', false);

descriptive_statistics = utils.calculate_descriptive_statistics(data, headers);

