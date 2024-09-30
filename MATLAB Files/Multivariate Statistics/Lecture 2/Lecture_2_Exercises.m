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

% Define column names based on the number of columns in your data
headers = arrayfun(@(x) sprintf('Dim %d', x), 1:size(data, 2), 'UniformOutput', false);

[mean, S, range] = utils.calculate_descriptive_statistics(data, headers);

% Display mean values
disp(array2table(mean, 'VariableNames', headers, 'RowNames', {'Mean'}));

fprintf('\nSariance Matrix (Diagonal):\n');
% Display Sariance matrix diagonal
disp(array2table(S, 'VariableNames', headers));

% Display range values
disp(array2table(range, 'VariableNames', headers, 'RowNames', {'Range'}));

% Mahalanobis
% d^2_i = (x_i - x_bar)^T * S^-1 * (x_i - x_bar)
mahalanobis_distances = [];
% data = data';

for j=1:size(data, 1)
    mahalanobis_distances = [mahalanobis_distances; (data(j, :) - mean) * inv(S) * (data(j, :) - mean)'];
end
% fprintf('Mahalanobis Distances Size: [%d %d]\n', size(mahalanobis_distances, 1), size(mahalanobis_distances, 2));
array2table(mahalanobis_distances, 'VariableNames', {'Mahalanobis Distances'})

figure(3)
qqplot(mahalanobis_distances, chi2rnd(ones(1, length(data(:, 1)))))

% The dataset does not appear linear like the chi squared distribution. It
% must be a different kind of distribution, like an exponential one.
%%

