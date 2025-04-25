clc;clear;close all
format compact
%--------------------------------------------------------------------------
% load data for ANOVA model with two factors of three levels each
%--------------------------------------------------------------------------
load('M:\Statistik\MULTISTAT\Lektion 06  UMLR\Course material\dataset_problem_6_3');
%--------------------------------------------------------------------------
% analyze ANOVA model with linear regression model:
% - build response vector Y
% - build design matrix Z including possible factor interactions
%--------------------------------------------------------------------------
Y = data(:,1);
Z = zeros(36,16);
% constant column 1
Z(:,1) = 1;
% column 2:4 for factor 1 (species)
Z(1:12,2) = 1;
Z(13:24,3) = 1;
Z(25:36,4) = 1;
% column 5:7 for factor 2 (time)
Z([1:4 13:16 25:28],5) = 1;
Z([5:8 17:20 29:32],6) = 1;
Z([9:12 21:24 33:36],7) = 1;
% column 8:16 for interactions beteween factor 1 and 2 
Z(1:4,8) = 1;
Z(5:8,9) = 1;
Z(9:12,10) = 1;
Z(13:16,11) = 1;
Z(17:20,12) = 1;
Z(21:24,13) = 1;
Z(25:28,14) = 1;
Z(29:32,15) = 1;
Z(33:36,16) = 1;
disp('-------------------------------------------------------------------------')
disp('Two-factor ANOVA with interactions as a special case of linear regression')
disp('-------------------------------------------------------------------------')
disp('ANOVA model:  mu_lk = mu + alfa_l + beta_k + gamma_lk,  l,k = 1,2,3 ')
disp('-------------------------------------------------------------------------')
disp('Multiple regression model:  Y = Zb + e')
disp('-------------------------------------------------------------------------')
Y_T = Y'
disp('-------------------------------------------------------------------------')
Z
disp('-------------------------------------------------------------------------')
rank_Z = rank(Z)
rank_ZT_Z = rank(Z'*Z)
det_ZT_Z = det(Z'*Z)
disp('-------------------------------------------------------------------------')
pseudo_inverse_ZT_Z = pinv(Z'*Z)
disp('-------------------------------------------------------------------------')
disp('Unbiased LS estimate of regression coefficients using pseudo-inverse')
disp('-------------------------------------------------------------------------')
b_hat = pseudo_inverse_ZT_Z*Z'*Y
disp('-------------------------------------------------------------------------')
disp('Extra Sums of squares based test for no interactions')
disp('H0: all gamma_lk = 0')
disp('-------------------------------------------------------------------------')
SSE_full_model = (Y - Z*b_hat)'*(Y - Z*b_hat)
n = 36
rank_Z = rank(Z)
df_full_model = n - rank_Z
sigma_hat_sqr = SSE_full_model/df_full_model
disp('-------------------------------------------------------------------------')
Z_no_interactions = Z(:,1:7);
b_hat_no_interactions = pinv(Z_no_interactions'*Z_no_interactions)*Z_no_interactions'*Y;
SSE_no_interactions_model = (Y - Z_no_interactions*b_hat_no_interactions)'*(Y - Z_no_interactions*b_hat_no_interactions)
rank_Z_no_interactions = rank(Z_no_interactions)
df_no_interactions_model = n - rank_Z_no_interactions
disp('-------------------------------------------------------------------------')
test_statistic = ((SSE_no_interactions_model - SSE_full_model)/(rank_Z - rank_Z_no_interactions))/sigma_hat_sqr
alpha = 0.05
critical_value = finv(1-alpha,rank_Z - rank_Z_no_interactions,df_full_model) 
rejection_of_H0 = test_statistic >= critical_value
p_value = 1-fcdf(test_statistic,rank_Z - rank_Z_no_interactions,df_full_model)
disp('-------------------------------------------------------------------------')
disp('Extra Sums of squares based test for no factor 1 effect (species)')
disp('H0: all alfa_l = 0')
disp('-------------------------------------------------------------------------')
Z_new = Z_no_interactions;
SSE_full_new_model = SSE_no_interactions_model;
n = 36
rank_Z_new = rank(Z_new)
df_full_new_model = n - rank_Z_new
sigma_hat_sqr_new = SSE_full_new_model/df_full_new_model
disp('-------------------------------------------------------------------------')
Z_no_factor1_effect = Z_new(:,[1 5:7]);
b_no_factor1_effect = pinv(Z_no_factor1_effect'*Z_no_factor1_effect)*Z_no_factor1_effect'*Y;
SSE_no_factor1_effect_model = (Y - Z_no_factor1_effect*b_no_factor1_effect)'*(Y - Z_no_factor1_effect*b_no_factor1_effect)
rank_Z_no_factor1_effect = rank(Z_no_factor1_effect)
disp('-------------------------------------------------------------------------')
test_statistic = ((SSE_no_factor1_effect_model - SSE_full_new_model)/(rank_Z_new - rank_Z_no_factor1_effect))/sigma_hat_sqr  % sigma_hat_sqr_new could be used
alpha = 0.05
critical_value = finv(1-alpha,rank_Z_new - rank_Z_no_factor1_effect,df_full_model) 
rejection_of_H0 = test_statistic >= critical_value
p_value = 1-fcdf(test_statistic,rank_Z_new - rank_Z_no_factor1_effect,df_full_model)
disp('-------------------------------------------------------------------------')
disp('Extra Sums of squares based test for no factor 2 effect (time)')
disp('H0: all beta_k = 0')
disp('-------------------------------------------------------------------------')
Z_new = Z_no_interactions;
SSE_full_new_model = SSE_no_interactions_model;
n = 36
rank_Z_new = rank(Z_new)
df_full_new_model = n - rank_Z_new
sigma_hat_sqr_new = SSE_full_new_model/df_full_new_model
disp('-------------------------------------------------------------------------')
Z_no_factor2_effect = Z_new(:,1:4);
b_no_factor2_effect = pinv(Z_no_factor2_effect'*Z_no_factor2_effect)*Z_no_factor2_effect'*Y;
SSE_no_factor2_effect_model = (Y - Z_no_factor2_effect*b_no_factor2_effect)'*(Y - Z_no_factor2_effect*b_no_factor2_effect)
rank_Z_no_factor2_effect = rank(Z_no_factor2_effect)
disp('-------------------------------------------------------------------------')
test_statistic = ((SSE_no_factor2_effect_model - SSE_full_new_model)/(rank_Z_new - rank_Z_no_factor2_effect))/sigma_hat_sqr  % sigma_hat_sqr_new could be used
alpha = 0.05
critical_value = finv(1-alpha,rank_Z_new - rank_Z_no_factor2_effect,df_full_model) 
rejection_of_H0 = test_statistic >= critical_value
p_value = 1-fcdf(test_statistic,rank_Z_new - rank_Z_no_factor2_effect,df_full_model)
disp('-------------------------------------------------------------------------')
disp('Comparison with ANOVA2 table')
disp('(anova2)')
disp('-------------------------------------------------------------------------')
X_anova2 = [Y(1:12) Y(13:24) Y(25:36)];
replications = 4;
[p_values_ANOVA2, ANOVA2_table] = anova2(X_anova2,replications,'off')
disp('-------------------------------------------------------------------------')
% %-------------------------------------------------------------------------------------------------------