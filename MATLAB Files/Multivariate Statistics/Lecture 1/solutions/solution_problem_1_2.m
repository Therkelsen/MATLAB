clc;clear;close all
format compact

%------------------------------------------------------------------------
% time data
%------------------------------------------------------------------------
X_time = load('C:\Users\cv\OneDrive - Syddansk Universitet\M-drive moved by SDU-IT\Statistik\MULTISTAT\Lektion 01  Multivariate Data\Course material\dataset_problem_1_2.dat');
disp('---------------------------------------------------------------------')
disp('Descriptive statistics estimators for female time-data')
disp('---------------------------------------------------------------------')
[n p] = size(X_time)
Xmean_time = mean(X_time,1)
S_time = cov(X_time)
R_time = corr(X_time)
figure('Name','Box plot for original female time-data')
boxplot(X_time)
figure('Name','Scatter plot for original female time-data')
plotmatrix(X_time,'g*')
generalized_variance_S_time = det(S_time)
generalized_variance_R_time = det(R_time)
prod_diag_S_generalized_variance_R_time = det(R_time)*prod(diag(S_time))

%------------------------------------------------------------------------
% speed data
%------------------------------------------------------------------------
disp('---------------------------------------------------------------------')
disp('Descriptive statistics estimators for female speed-data')
disp('---------------------------------------------------------------------')
X_speed(:,1) = 100./X_time(:,1);
X_speed(:,2) = 200./X_time(:,2);
X_speed(:,3) = 400./X_time(:,3);
X_speed(:,4) = 800./(60*X_time(:,4));
X_speed(:,5) = 1500./(60*X_time(:,5));
X_speed(:,6) = 3000./(60*X_time(:,6));
X_speed(:,7) = 42195./(60*X_time(:,7));
Xmean_speed = mean(X_speed,1)
S_speed = cov(X_speed)
R_speed = corr(X_speed)
figure('Name','Box plot for original female speed-data')
boxplot(X_speed)
figure('Name','Scatter plot for original female speed-data')
plotmatrix(X_speed,'r*')
generalized_variance_S_speed = det(S_speed)
generalized_variance_R_speed = det(R_speed)
prod_diag_S_generalized_variance_R_speed = det(R_speed)*prod(diag(S_speed))
disp('---------------------------------------------------------------------')

%--------------------------------------------------------------------------