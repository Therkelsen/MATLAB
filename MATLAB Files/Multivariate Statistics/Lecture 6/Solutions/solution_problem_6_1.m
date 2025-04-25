clc;clear;close all
format compact
%--------------------------------------------------------------------------
% load data
%--------------------------------------------------------------------------
% manual calculation
%--------------------------------------------------------------------------
load('C:\Users\cv\OneDrive - Syddansk Universitet\M-drive moved by SDU-IT\Statistik\MULTISTAT\Lektion 06  UMLR\Course material\dataset_problem_6_1');
disp('-------------------------------------------------------------------')
disp('Multiple regression model:  Y = b0 + b1*z1 + b2*z2 + e')
disp('-------------------------------------------------------------------')
n = length(Y)
r = 2
Z = [ones(n,1) Z];
disp('-------------------------------------------------------------------')
disp('Unbiased Least Squares estimate of regression coefficients')
disp('-------------------------------------------------------------------')
b_hat = inv(Z'*Z)*Z'*Y;
b0_hat = b_hat(1)
b1_hat = b_hat(2)
b2_hat = b_hat(3)
disp('-------------------------------------------------------------------')
disp('Unbiased estimate of Var[e]')
disp('-------------------------------------------------------------------')
Y_hat = Z*b_hat;
e_hat = Y - Y_hat;
SSE = e_hat'*e_hat;
sigma_hat_sqr = SSE/(n-(r+1))
sigma_hat = sqrt(sigma_hat_sqr)
disp('-------------------------------------------------------------------')
disp('SS decomposition: SST = SSR + SSE')
disp('-------------------------------------------------------------------')
SST = (Y-mean(Y))'*(Y-mean(Y))
SSR = (Y_hat-mean(Y))'*(Y_hat-mean(Y))
SSE
disp('-------------------------------------------------------------------')
disp('Coefficient of determination, R_sqr, and R_adj_sqr')
disp('-------------------------------------------------------------------')
R_sqr = 1 - SSE/SST
R_adj_sqr = 1 - (SSE/(n-(r+1)))/(SST/(n-1))
disp('-------------------------------------------------------------------')
disp('Estimated covariance matrix for b_hat')
disp('-------------------------------------------------------------------')
SIGMA_hat_b_hat = sigma_hat_sqr*inv(Z'*Z)
sigma_hat_b0_hat = sqrt(SIGMA_hat_b_hat(1,1))
sigma_hat_b1_hat = sqrt(SIGMA_hat_b_hat(2,2))
sigma_hat_b2_hat = sqrt(SIGMA_hat_b_hat(3,3))
disp('-------------------------------------------------------------------')
disp('95% simultaneous confidence intervals for b')
disp('-------------------------------------------------------------------')
alpha = 0.05;
sim_CI_b0 = b0_hat + [-1 1]*sqrt((r+1)*finv(1-alpha,r+1,n-(r+1)))*sigma_hat_b0_hat
sim_CI_b1 = b1_hat + [-1 1]*sqrt((r+1)*finv(1-alpha,r+1,n-(r+1)))*sigma_hat_b1_hat
sim_CI_b2 = b2_hat + [-1 1]*sqrt((r+1)*finv(1-alpha,r+1,n-(r+1)))*sigma_hat_b2_hat
disp('-------------------------------------------------------------------')
disp('95% Bonferroni confidence intervals for b')
disp('-------------------------------------------------------------------')
bonf_CI_b0 = b0_hat + [-1 1]*tinv(1-alpha/(2*(r+1)),n-(r+1))*sigma_hat_b0_hat
bonf_CI_b1 = b1_hat + [-1 1]*tinv(1-alpha/(2*(r+1)),n-(r+1))*sigma_hat_b1_hat
bonf_CI_b2 = b2_hat + [-1 1]*tinv(1-alpha/(2*(r+1)),n-(r+1))*sigma_hat_b2_hat
disp('-------------------------------------------------------------------')
disp('Extra Sum of Squares based test for significant dependency of Y on (z1,z2)')
disp('(significant regression model)')
disp('H0: b1 = b2 = 0, b0 arbitrary')
disp('-------------------------------------------------------------------')
q = 0
SSE_full_model = (Y - Z*b_hat)'*(Y - Z*b_hat)
Z_only_intercept = Z(:,1:(q+1));
b_hat_only_intercept = inv(Z_only_intercept'*Z_only_intercept)*Z_only_intercept'*Y
SSE_only_intercept_model = (Y - Z_only_intercept*b_hat_only_intercept)'*(Y - Z_only_intercept*b_hat_only_intercept)
test_statistic = ((SSE_only_intercept_model - SSE_full_model)/(r-q))/sigma_hat_sqr
critical_value = finv(1-alpha,r-q,n-(r+1)) 
rejection_of_H0 = test_statistic >= critical_value
p_value = 1-fcdf(test_statistic,r-q,n-(r+1))
%--------------------------------------------------------------------------
% forkert metode :  (ny estimation af sigma_sqr)  samme konklusion dog
sigma_hat_sqr_only_intercept_model = SSE_only_intercept_model/(n-(q+1))
test_statistic_alt = ((SSE_only_intercept_model - SSE_full_model)/(r-q))/sigma_hat_sqr_only_intercept_model
critical_value_alt = finv(1-alpha,r-q,n-(q+1)) 
rejection_of_H0_alt = test_statistic_alt >= critical_value_alt
p_value_alt = 1-fcdf(test_statistic_alt,r-q,n-(q+1))
disp('-------------------------------------------------------------------')
disp('ANOVA table based test for significant dependency of Y on (z1,z2)')
disp('(significant regression model)')
disp('H0: b1 = b2 = 0, b0 arbitrary')
disp('-------------------------------------------------------------------')
SSR
df_SSR = r
SSE
df_SSE = n-(r+1)
test_statistic = (SSR/df_SSR)/(SSE/df_SSE)
critical_value = finv(1-alpha,df_SSR,df_SSE) 
rejection_of_H0 = test_statistic >= critical_value
p_value = 1-fcdf(test_statistic,df_SSR,df_SSE)
disp('-------------------------------------------------------------------')
disp('Marginal test for  H0: b0 = 0')
disp('-------------------------------------------------------------------')
b0_hat
sigma_hat_b0_hat
test_statistic = b0_hat/sigma_hat_b0_hat
critical_value = tinv(1-alpha/2,n-(r+1))
rejection_of_H0 = test_statistic >= critical_value
p_value = 2*(1-tcdf(test_statistic,n-(r+1)))
disp('-------------------------------------------------------------------')
disp('Marginal test for  H0: b1 = 0')
disp('-------------------------------------------------------------------')
b1_hat
sigma_hat_b1_hat
test_statistic = b1_hat/sigma_hat_b1_hat
critical_value = tinv(1-alpha/2,n-(r+1))
rejection_of_H0 = test_statistic >= critical_value
p_value = 2*(1-tcdf(test_statistic,n-(r+1)))
disp('-------------------------------------------------------------------')
disp('Marginal test for  H0: b2 = 0')
disp('-------------------------------------------------------------------')
b2_hat
sigma_hat_b2_hat
test_statistic = b2_hat/sigma_hat_b2_hat
critical_value = tinv(1-alpha/2,n-(r+1))
rejection_of_H0 = test_statistic >= critical_value
p_value = 2*(1-tcdf(test_statistic,n-(r+1)))
disp('-------------------------------------------------------------------')
disp('95% prediction interval for E[Y(z0)] (mean of new responses at z=z0)')
disp('-------------------------------------------------------------------')
z0_T = [1 50.5 970]
z0 = z0_T';
Y_z0_hat = z0'*b_hat
pred_CI_mean_Y_z0 = Y_z0_hat + [-1 1]*sqrt(z0'*inv(Z'*Z)*z0)*tinv(1-alpha/2,n-(r+1))*sigma_hat
disp('-------------------------------------------------------------------')
disp('95% prediction interval for Y(z0) (single new response at z=z0)')
disp('-------------------------------------------------------------------')
z0_T = [1 50.5 970]
z0 = z0_T';
Y_z0_hat = z0'*b_hat
pred_CI_Y_z0 = Y_z0_hat + [-1 1]*sqrt(1+z0'*inv(Z'*Z)*z0)*tinv(1-alpha/2,n-(r+1))*sigma_hat
disp('-------------------------------------------------------------------')
disp('Model check: analysis of residuals (estimated model errors)')
disp('SEE FIGURE 1')
disp('-------------------------------------------------------------------')
e_hat = Y - Y_hat;
figure(1)
subplot(7,3,1:3)
axis off
text(0.04,0.4,'Analysis of residuals (estimated model errors),  e_j_,_h_a_t = Y_j - Y_j_,_h_a_t = Y_j - Z*b_h_a_t,  j = 1,2,...,n','Fontsize',20)
subplot(7,3,[7 10])
j = 1:n;
plot(j,e_hat,'r+'),line([1 n],[0 0],'color','k'),xlim([1 n]),title('e_j_,_h_a_t = f(j)','Fontsize',14)
subplot(7,3,[8 11])
plot(Y,e_hat,'k+'),line([min(Y) max(Y)],[0 0],'color','k'),xlim([min(Y) max(Y)]),title('e_j_,_h_a_t = f(Y_j)','Fontsize',14)
subplot(7,3,[16 19])
Z1 = Z(:,2);
plot(Z1,e_hat,'g+'),line([min(Z1) max(Z1)],[0 0],'color','k'),xlim([min(Z1) max(Z1)]),title('e_j_,_h_a_t = f(z_j_1)','Fontsize',14)
subplot(7,3,[17 20])
Z2 = Z(:,3);
plot(Z2,e_hat,'g+'),line([min(Z2) max(Z2)],[0 0],'color','k'),xlim([min(Z2) max(Z2)]),title('e_j_,_h_a_t = f(z_j_2)','Fontsize',14)
subplot(7,3,[9 12])
hist(e_hat),title('histogram of e_j_,_h_a_t','Fontsize',14)
subplot(7,3,[18 21])
normplot(e_hat),title('normal probability plot (QQ-plot)','Fontsize',14)
%--------------------------------------------------------------------------
% calculations using built-in functions
%--------------------------------------------------------------------------
disp(' ')
disp('-----------------------------------------------------------------------------------------------------------')
disp('-----------------------------------------------------------------------------------------------------------')
disp(' ')
disp('-------------------------------------------------------------------')
disp('Using built-in function:   ')
disp('LinearModel.fit --->  create model')
disp('LinearModel.fit --->  estimate of regression coefficients')
disp('LinearModel.fit --->  marginal tests for coefficients')
disp('LinearModel.fit --->  R-squared and adjusted R-squared')
disp('LinearModel.fit --->  test for significant regression model')
disp('-------------------------------------------------------------------')
fitted_model = LinearModel.fit(Z(:,2:3),Y)
disp('-------------------------------------------------------------------')
disp('Using built-in function:   ')
disp('anova --->  marginal tests for coefficients (not intercept)')
disp('-------------------------------------------------------------------')
anova_table = anova(fitted_model)
disp('-------------------------------------------------------------------')
disp('Using built-in methods in class LinearModel:   ')
disp('coefCI --->  95% marginal confidence intervals for coefficients')
disp('-------------------------------------------------------------------')
marg_CI_b = coefCI(fitted_model)
disp('-------------------------------------------------------------------')
disp('Using built-in methods in class LinearModel:   ')
disp('predict --->  95% confidence intervals for new observations')
disp('-------------------------------------------------------------------')
z0 = [50.5 970]
[Y_z0_hat, pred_CI_mean_Y_z0]= predict(fitted_model,z0,'Prediction','Curve')
[Y_z0_hat, pred_CI_Y_z0]= predict(fitted_model,z0,'Prediction','Observation');
pred_CI_Y_z0
disp('-------------------------------------------------------------------')
disp('Using built-in methods in class LinearModel:   ')
disp('plotDiagnostics --->  SEE FIGURE 2')
disp('plotResiduals   --->  SEE FIGURE 2')
disp('-------------------------------------------------------------------')
figure(2)
subplot(2,2,1)
plotDiagnostics(fitted_model)
subplot(2,2,2)
plotResiduals(fitted_model)
subplot(2,2,3)
plotResiduals(fitted_model,'probability')
subplot(2,2,4)
plotResiduals(fitted_model,'fitted')
%-------------------------------------------------------------------------------------------------------