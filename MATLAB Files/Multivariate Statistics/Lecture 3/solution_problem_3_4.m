clc;clear;close all
format compact
%--------------------------------------------------------------------------
disp('-----------------------------------------------------')
disp('Likelihood Ratio Test for Poisson distribution')
disp('H0: lambda = lambda0')
lambda0 = 10
n = 100
lambda = lambda0
%lambda = 8
x_i = poissrnd(lambda,1,n);
x_mean = mean(x_i)
test_statistic = 2*n*(lambda0 - x_mean - x_mean*log(lambda0/x_mean))
alfa = 0.05
df = 1
critical_value = chi2inv(1-alfa,df)
rejction_of_H0 = test_statistic >= critical_value
p_value = 1 - chi2cdf(test_statistic,df)
%--------------------------------------------------------------------------
disp('-----------------------------------------------------')
disp('Likelihood Ratio Test for exponential distribution')
disp('H0: lambda = lambda0')
lambda0 = 10
n = 100
lambda = lambda0
%lambda = 8
t_i = exprnd(1/lambda,1,n);
t_mean = mean(t_i)
test_statistic = 2*n*(lambda0*t_mean - log(lambda0*t_mean) - 1)
alfa = 0.05
df = 1
critical_value = chi2inv(1-alfa,df)
rejction_of_H0 = test_statistic >= critical_value
p_value = 1 - chi2cdf(test_statistic,df)
disp('-----------------------------------------------------')
%--------------------------------------------------------------------------
