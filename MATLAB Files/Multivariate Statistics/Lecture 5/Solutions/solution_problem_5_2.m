clc;clear;close all
format compact
%--------------------------------------------------------------------------
% load data
%--------------------------------------------------------------------------
load('C:\Users\cv\OneDrive - Syddansk Universitet\M-drive moved by SDU-IT\Statistik\MULTISTAT\Lektion 05  MANOVA\Course material\dataset_problem_5_2');
disp('---------------------------------------------------------------------------')
disp('MANOVA 2')
disp('---------------------------------------------------------------------------')
%--------------------------------------------------------------------------
% data parameters
%--------------------------------------------------------------------------
disp('data parameters')
disp('---------------------------------------------------------------------------')
p = 3
g = 2
b = 3
n = 2
N = g*b*n
%--------------------------------------------------------------------------
% sample means 
%--------------------------------------------------------------------------
disp('---------------------------------------------------------------------------')
disp('sample means for factor 1+2 combinations')
disp('---------------------------------------------------------------------------')
X11 = X(1:2,3:5);
mu_hat_11 = mean(X11)
X21 = X(3:4,3:5);
mu_hat_21 = mean(X21)
X12 = X(5:6,3:5);
mu_hat_12 = mean(X12)
X22 = X(7:8,3:5);
mu_hat_22 = mean(X22)
X13 = X(9:10,3:5);
mu_hat_13 = mean(X13)
X23 = X(11:12,3:5);
mu_hat_23 = mean(X23)
disp('---------------------------------------------------------------------------')
disp('sample means for factor 1')
disp('---------------------------------------------------------------------------')
mu_hat_1x = mean([mu_hat_11; mu_hat_12; mu_hat_13])
mu_hat_2x = mean([mu_hat_21; mu_hat_22; mu_hat_23])
disp('---------------------------------------------------------------------------')
disp('sample means for factor 2')
disp('---------------------------------------------------------------------------')
mu_hat_x1 = mean([mu_hat_11; mu_hat_21])
mu_hat_x2 = mean([mu_hat_12; mu_hat_22])
mu_hat_x3 = mean([mu_hat_13; mu_hat_23])
disp('---------------------------------------------------------------------------')
disp('overall sample mean')
disp('---------------------------------------------------------------------------')
mu_hat = mean(X(:,3:5))
%--------------------------------------------------------------------------
% SS decomposition
%--------------------------------------------------------------------------
disp('---------------------------------------------------------------------------')
disp('SS decomposition')
disp('---------------------------------------------------------------------------')
SSB_factor1_1 = b*n*(mu_hat_1x - mu_hat)'*(mu_hat_1x - mu_hat);
SSB_factor1_2 = b*n*(mu_hat_2x - mu_hat)'*(mu_hat_2x - mu_hat);
SSB_factor1 = SSB_factor1_1 + SSB_factor1_2
df_factor1 = g-1
SB_factor1 = SSB_factor1/df_factor1
disp('---------------------------------------------------------------------------')
SSB_factor2_1 = g*n*(mu_hat_x1 - mu_hat)'*(mu_hat_x1 - mu_hat);
SSB_factor2_2 = g*n*(mu_hat_x2 - mu_hat)'*(mu_hat_x2 - mu_hat);
SSB_factor2_3 = g*n*(mu_hat_x3 - mu_hat)'*(mu_hat_x3 - mu_hat);
SSB_factor2 = SSB_factor2_1 + SSB_factor2_2 + SSB_factor2_3
df_factor2 = b-1
SB_factor2 = SSB_factor2/df_factor2
disp('---------------------------------------------------------------------------')
SS_interact_11 = n*(mu_hat_11 - mu_hat_1x - mu_hat_x1 + mu_hat)'*(mu_hat_11 - mu_hat_1x - mu_hat_x1 + mu_hat);
SS_interact_21 = n*(mu_hat_21 - mu_hat_2x - mu_hat_x1 + mu_hat)'*(mu_hat_21 - mu_hat_2x - mu_hat_x1 + mu_hat);
SS_interact_12 = n*(mu_hat_12 - mu_hat_1x - mu_hat_x2 + mu_hat)'*(mu_hat_12 - mu_hat_1x - mu_hat_x2 + mu_hat);
SS_interact_22 = n*(mu_hat_22 - mu_hat_2x - mu_hat_x2 + mu_hat)'*(mu_hat_22 - mu_hat_2x - mu_hat_x2 + mu_hat);
SS_interact_13 = n*(mu_hat_13 - mu_hat_1x - mu_hat_x3 + mu_hat)'*(mu_hat_13 - mu_hat_1x - mu_hat_x3 + mu_hat);
SS_interact_23 = n*(mu_hat_23 - mu_hat_2x - mu_hat_x3 + mu_hat)'*(mu_hat_23 - mu_hat_2x - mu_hat_x3 + mu_hat);
SS_interact = SS_interact_11 + SS_interact_21 + SS_interact_12 + SS_interact_22 + SS_interact_13 + SS_interact_23
df_interact = (g-1)*(b-1)
S_interact = SS_interact/df_interact
disp('---------------------------------------------------------------------------')
SS_total = (X(:,3:5) - ones(N,1)*mu_hat)'*(X(:,3:5) - ones(N,1)*mu_hat);
SSW = SS_total - SSB_factor1 - SSB_factor2 - SS_interact
df_W = g*b*(n-1)
SW = SSW/df_W
disp('---------------------------------------------------------------------------')
SS_total
df_total = g*b*n-1
S_total = SS_total/df_total
%--------------------------------------------------------------------------
% MANOVA2 test for no systematic interaction (all gamma_lk = 0)
%--------------------------------------------------------------------------
disp('--------------------------------------------------------------------------')
disp('MANOVA2 test for no systematic interaction (all gamma_lk = 0)')
disp('--------------------------------------------------------------------------')
LAMBDA = det(SSW)/det(SS_interact + SSW)
test_statistic = -(g*b*(n-1) -((p+1)-(g-1)*(b-1))/2)*log(LAMBDA)
alpha = 0.01
critical_value = chi2inv(1-alpha,p*(g-1)*(b-1))
rejction_of_H0 = test_statistic >= critical_value
p_value = 1-chi2cdf(test_statistic,p*(g-1)*(b-1))
%--------------------------------------------------------------------------
% MANOVA2 test for no factor 1 effect (all tau_1 = 0)
%--------------------------------------------------------------------------
disp('--------------------------------------------------------------------------')
disp('MANOVA2 test for no factor 1 effect (all tau_1 = 0)')
disp('--------------------------------------------------------------------------')
LAMBDA = det(SSW)/det(SSB_factor1 + SSW)
test_statistic = -(g*b*(n-1) -((p+1)-(g-1))/2)*log(LAMBDA)
alpha = 0.01
critical_value = chi2inv(1-alpha,p*(g-1))
rejction_of_H0 = test_statistic >= critical_value
p_value = 1-chi2cdf(test_statistic,p*(g-1))
%--------------------------------------------------------------------------
% MANOVA2 test for no factor 2 effect (all beta_k = 0)
%--------------------------------------------------------------------------
disp('--------------------------------------------------------------------------')
disp('MANOVA2 test for no factor 2 effect (all beta_k = 0)')
disp('--------------------------------------------------------------------------')
LAMBDA = det(SSW)/det(SSB_factor2 + SSW)
test_statistic = -(g*b*(n-1) -((p+1)-(b-1))/2)*log(LAMBDA)
alpha = 0.01
critical_value = chi2inv(1-alpha,p*(b-1))
rejction_of_H0 = test_statistic >= critical_value
p_value = 1-chi2cdf(test_statistic,p*(b-1))
%--------------------------------------------------------------------------
% CONCLUSION:  X_lkr ~ N_3(mu_k,SIGMA)
%--------------------------------------------------------------------------
disp('--------------------------------------------------------------------------')
disp('CONCLUSION:  X_lkr ~ N_3(mu_k,SIGMA) with parameter estimates:')
disp('--------------------------------------------------------------------------')
mu_1_hat = mu_hat_x1
mu_2_hat = mu_hat_x2
mu_3_hat = mu_hat_x3
SIGMA_hat = (SS_total - SSB_factor2)/(df_total - df_factor2)
disp('--------------------------------------------------------------------------')
%--------------------------------------------------------------------------
