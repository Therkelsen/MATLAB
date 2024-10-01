clc; close all;
format compact

% Exercise 1.1
rho = 28;
sigma = 10;
beta = 8/3;

% Exercise 1.2
la_a = out.lorenz;
plot3(la_a(:, 1), la_a(:, 2), la_a(:, 3))