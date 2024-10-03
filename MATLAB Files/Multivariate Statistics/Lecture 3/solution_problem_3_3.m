clc;clear;close all
format compact
%--------------------------------------------------------------------------
disp('-----------------------------------------------------')
disp('Likelihood Ratio Test for multivariate normal distribution')
disp('H0: SIGMA = SIGMA_0')
p = 2
mu = [2 4]
SIGMA_0 = [1 1.5; 1.5 4]
n = 100
SIGMA = SIGMA_0
%SIGMA = [2 1.5; 1.5 4]
x_i = mvnrnd(mu,SIGMA,n);
S = cov(x_i)
test_statistic = n*trace(inv(SIGMA_0)*S) - n*log(det(inv(SIGMA_0)*S)) - n*p
alfa = 0.05
df = p*(p+1)/2
critical_value = chi2inv(1-alfa,df)
rejction_of_H0 = test_statistic >= critical_value
p_value = 1 - chi2cdf(test_statistic,df)
%--------------------------------------------------------------------------
