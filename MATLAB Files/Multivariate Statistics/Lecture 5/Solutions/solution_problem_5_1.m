clc;clear;close all
format compact
%--------------------------------------------------------------------------
% load data
%--------------------------------------------------------------------------
load('C:\Users\cv\OneDrive - Syddansk Universitet\M-drive moved by SDU-IT\Statistik\MULTISTAT\Lektion 05  MANOVA\Course material\dataset_problem_5_1');
disp('-------------------------------------------------------------------')
disp('Bartlett test,  H0: SIGMA_time1 = SIGMA_time2 = SIGMA_time3')
disp('-------------------------------------------------------------------')
%--------------------------------------------------------------------------
% descriptive statistics for Bartlett test
%--------------------------------------------------------------------------
p = 4
g = 3
n_time1 = 30
n_time2 = 30
n_time3 = 30
n = n_time1 + n_time2 + n_time3
Sigma_hat_time1 = cov(X_time1)
Sigma_hat_time2 = cov(X_time2)
Sigma_hat_time3 = cov(X_time3)
Sigma_hat_pooled = ((n_time1-1)*Sigma_hat_time1 + (n_time2-1)*Sigma_hat_time2 + (n_time3-1)*Sigma_hat_time3)/(n-g)
%--------------------------------------------------------------------------
%  Bartlett test (Box's M-test) for equal covariance matrices
%--------------------------------------------------------------------------
T = (n-g)*log(det(Sigma_hat_pooled)) - (n_time1-1)*log(det(Sigma_hat_time1)) - (n_time2-1)*log(det(Sigma_hat_time2)) - (n_time3-1)*log(det(Sigma_hat_time3))
correction_factor = 1 - ((2*p^2+3*p-1)/(6*(p+1)*(g-1)))*(1/(n_time1-1)+1/(n_time2-1)+1/(n_time3-1)-1/(n-g))
test_statistic = correction_factor*T
alpha = 0.05
df = p*(p+1)*(g-1)/2
critical_value = chi2inv(1-alpha,df)
rejection_of_H0 = test_statistic >= critical_value
p_value = 1-chi2cdf(test_statistic,df)
%[ndim, prob] = barttest([X_time1 X_time2 X_time3],alpha)
%------------------------------------------------------------------------
% Model check for TIME 1 data
%------------------------------------------------------------------------
figure
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for TIME 1 observations possibly being from multivariate normal distribution (~ MVN with p = 4)','Fontsize',20)
subplot(5,2,3:2:9)
plotmatrix(X_time1,'b+')
title('Scattermatrix of TIME 1 observations','Fontsize',16)
subplot(5,2,4)
axis off
text(0,1,'Sample Mahalanobis distances, i=1..n_t_i_m_e_1 :','Fontsize',18)
text(0,0.7,'d_i^2 = (x_i - x)^T S_t_i_m_e_1^-^1 (x_i - x) ~ {\chi}_4^2   (approximately)','Fontsize',16)
text(0.158,0.91,'_','Fontsize',16),warning('off')
text(0.4,0.91,'_','Fontsize',16),warning('off')
subplot(5,2,6:2:10)
mu_hat_time1 = mean(X_time1);
d_i_sqr = zeros(1,n_time1);
for i = 1:n_time1
    d_i_sqr(i) = (X_time1(i,:)-mu_hat_time1)*inv(Sigma_hat_time1)*(X_time1(i,:)-mu_hat_time1)';
end
df = p;
z_i = chi2rnd(df,1,n_time1);
% z2 = chi2rnd(20,1,n);
qqplot(d_i_sqr,z_i),grid,xlabel('quantiles for d_i^2','Fontsize',16),ylabel('quantiles for {\chi}_4^2 distribution','Fontsize',16),...
    title('qq-plot for d_i^2 versus {\chi}_4^2 distribution','Fontsize',16)
%------------------------------------------------------------------------
% Model check for TIME 2 data
%------------------------------------------------------------------------
figure
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for TIME 2 observations possibly being from multivariate normal distribution (~ MVN with p = 4)','Fontsize',20)
subplot(5,2,3:2:9)
plotmatrix(X_time2,'gx')
title('Scattermatrix of TIME 2 observations','Fontsize',16)
subplot(5,2,4)
axis off
text(0,1,'Sample Mahalanobis distances, i=1..n_t_i_m_e_2 :','Fontsize',18)
text(0,0.7,'d_i^2 = (x_i - x)^T S_t_i_m_e_2^-^1 (x_i - x) ~ {\chi}_4^2   (approximately)','Fontsize',16)
text(0.158,0.91,'_','Fontsize',16),warning('off')
text(0.4,0.91,'_','Fontsize',16),warning('off')
subplot(5,2,6:2:10)
mu_hat_time2 = mean(X_time2);
d_i_sqr = zeros(1,n_time2);
for i = 1:n_time2
    d_i_sqr(i) = (X_time2(i,:)-mu_hat_time2)*inv(Sigma_hat_time2)*(X_time2(i,:)-mu_hat_time2)';
end
df = p;
z_i = chi2rnd(df,1,n_time2);
% z2 = chi2rnd(20,1,n);
qqplot(d_i_sqr,z_i),grid,xlabel('quantiles for d_i^2','Fontsize',16),ylabel('quantiles for {\chi}_4^2 distribution','Fontsize',16),...
    title('qq-plot for d_i^2 versus {\chi}_4^2 distribution','Fontsize',16)
%------------------------------------------------------------------------
% Model check for TIME 3 data
%------------------------------------------------------------------------
figure
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for TIME 3 observations possibly being from multivariate normal distribution (~ MVN with p = 4)','Fontsize',20)
subplot(5,2,3:2:9)
plotmatrix(X_time3,'ro')
title('Scattermatrix of TIME 3 observations','Fontsize',16)
subplot(5,2,4)
axis off
text(0,1,'Sample Mahalanobis distances, i=1..n_t_i_m_e_3 :','Fontsize',18)
text(0,0.7,'d_i^2 = (x_i - x)^T S_t_i_m_e_3^-^1 (x_i - x) ~ {\chi}_4^2   (approximately)','Fontsize',16)
text(0.158,0.91,'_','Fontsize',16),warning('off')
text(0.4,0.91,'_','Fontsize',16),warning('off')
subplot(5,2,6:2:10)
mu_hat_time3 = mean(X_time3);
d_i_sqr = zeros(1,n_time3);
for i = 1:n_time3
    d_i_sqr(i) = (X_time3(i,:)-mu_hat_time3)*inv(Sigma_hat_time3)*(X_time3(i,:)-mu_hat_time3)';
end
df = p;
z_i = chi2rnd(df,1,n_time3);
% z2 = chi2rnd(20,1,n);
qqplot(d_i_sqr,z_i),grid,xlabel('quantiles for d_i^2','Fontsize',16),ylabel('quantiles for {\chi}_4^2 distribution','Fontsize',16),...
    title('qq-plot for d_i^2 versus {\chi}_4^2 distribution','Fontsize',16)
%--------------------------------------------------------------------------
% Plot of observations from all g = 3 groups together in one scatterplot
%--------------------------------------------------------------------------
X_time = [X_time1; X_time2; X_time3];
time_index = [1*ones(n_time1,1); 2*ones(n_time2,1); 3*ones(n_time3,1)];
figure
gplotmatrix(X_time,[],time_index,[],'+xo'),title('observations from TIME 1, 2, 3  plotted together in one scatterplot','Fontsize',18)
%--------------------------------------------------------------------------
% Exact MANOVA1 test of factor effect (manual)
%--------------------------------------------------------------------------
disp('-------------------------------------------------------------------')
disp('Exact MANOVA1 test,  H0: mu_time1 = mu_time2 = mu_time3')
disp('-------------------------------------------------------------------')
mu_hat = mean([mu_hat_time1; mu_hat_time2; mu_hat_time3]);
SSB_time1 = n_time1*(mu_hat_time1 - mu_hat)'*(mu_hat_time1 - mu_hat);
SSB_time2 = n_time2*(mu_hat_time2 - mu_hat)'*(mu_hat_time2 - mu_hat);
SSB_time3 = n_time3*(mu_hat_time3 - mu_hat)'*(mu_hat_time3 - mu_hat);
SSB = SSB_time1 + SSB_time2 + SSB_time3
df_B = g-1
SB = SSB/df_B
SSW_time1 = (X_time1 - ones(n_time1,1)*mu_hat_time1)'*(X_time1 - ones(n_time1,1)*mu_hat_time1);
SSW_time2 = (X_time2 - ones(n_time2,1)*mu_hat_time2)'*(X_time2 - ones(n_time2,1)*mu_hat_time2);
SSW_time3 = (X_time3 - ones(n_time3,1)*mu_hat_time3)'*(X_time3 - ones(n_time3,1)*mu_hat_time3);
SSW = SSW_time1 + SSW_time2 + SSW_time3
df_W = n-g
SW = SSW/df_W
LAMBDA = det(SSW)/det(SSB + SSW);
test_statistic = ((n-p-2)/p)*(1-sqrt(LAMBDA))/sqrt(LAMBDA)
alpha = 0.05
critical_value = finv(1-alpha,2*p,2*(n-p-2))
rejection_of_H0 = test_statistic >= critical_value
p_value = 1-fcdf(test_statistic,2*p,2*(n-p-2))
%--------------------------------------------------------------------------
% Approximate large-sample MANOVA1 test of factor effect (manual)
%--------------------------------------------------------------------------
disp('--------------------------------------------------------------------------')
disp('Approximate large-sample MANOVA1 test,  H0: mu_time1 = mu_time2 = mu_time3')
disp('--------------------------------------------------------------------------')
test_statistic = -(n-1-(p+g)/2)*log(LAMBDA)
alpha = 0.05
critical_value = chi2inv(1-alpha,p*(g-1))
rejection_of_H0 = test_statistic >= critical_value
p_value = 1-chi2cdf(test_statistic,p*(g-1))
%--------------------------------------------------------------------------
% MANOVA1 test of factor effect with "manova1"
%--------------------------------------------------------------------------
disp('--------------------------------------------------------------------------------')
disp('MANOVA1 test of factor effect with "manova1",  H0: mu_time1 = mu_time2 = mu_time3')
disp('--------------------------------------------------------------------------------')
[dim,p_values,stats] = manova1(X_time,time_index);
SSB = stats.B
SSW = stats.W
df_B = stats.dfB
df_W = stats.dfW
chisq = stats.chisq(1)
p_value = p_values(1)
%--------------------------------------------------------------------------
% Parameter estimates for accept of no factor effect
%--------------------------------------------------------------------------
disp('--------------------------------------------------------------------------')
disp('Parameter estimates for accept of no factor effect')
disp('--------------------------------------------------------------------------')
mu_hat
SIGMA_hat = (SSB+SSW)/(n-1) 
%--------------------------------------------------------------------------
% Parameter estimates for rejection of no factor effect
%--------------------------------------------------------------------------
disp('--------------------------------------------------------------------------')
disp('Parameter estimates for rejection of no factor effect')
disp('--------------------------------------------------------------------------')
mu_hat_time1
mu_hat_time2
mu_hat_time3
SIGMA_hat = SW
%--------------------------------------------------------------------------
% Pairwise, simultaneous Bonferroni confidence intervals for difference in means
%--------------------------------------------------------------------------
disp('------------------------------------------------------------------------------')
disp('Pairwise, simultaneous Bonferroni confidence intervals for difference in means')
disp('------------------------------------------------------------------------------')
alpha = 0.05
disp('-------------')
disp('Variable X1')
disp('-------------')
bonf_CI_mu_diff_12_X1 = mu_hat_time1(1)-mu_hat_time2(1) + [-1 1]*tinv(1-alpha/(p*g*(g-1)),n-g)*sqrt(SIGMA_hat(1,1)*(1/n_time1+1/n_time2))
bonf_CI_mu_diff_13_X1 = mu_hat_time1(1)-mu_hat_time3(1) + [-1 1]*tinv(1-alpha/(p*g*(g-1)),n-g)*sqrt(SIGMA_hat(1,1)*(1/n_time1+1/n_time3))
bonf_CI_mu_diff_23_X1 = mu_hat_time2(1)-mu_hat_time3(1) + [-1 1]*tinv(1-alpha/(p*g*(g-1)),n-g)*sqrt(SIGMA_hat(1,1)*(1/n_time2+1/n_time3))
disp('-------------')
disp('Variable X2')
disp('-------------')
bonf_CI_mu_diff_12_X2 = mu_hat_time1(2)-mu_hat_time2(2) + [-1 1]*tinv(1-alpha/(p*g*(g-1)),n-g)*sqrt(SIGMA_hat(2,2)*(1/n_time1+1/n_time2))
bonf_CI_mu_diff_13_X2 = mu_hat_time1(2)-mu_hat_time3(2) + [-1 1]*tinv(1-alpha/(p*g*(g-1)),n-g)*sqrt(SIGMA_hat(2,2)*(1/n_time1+1/n_time3))
bonf_CI_mu_diff_23_X2 = mu_hat_time2(2)-mu_hat_time3(2) + [-1 1]*tinv(1-alpha/(p*g*(g-1)),n-g)*sqrt(SIGMA_hat(2,2)*(1/n_time2+1/n_time3))
disp('-------------')
disp('Variable X3')
disp('-------------')
bonf_CI_mu_diff_12_X3 = mu_hat_time1(3)-mu_hat_time2(3) + [-1 1]*tinv(1-alpha/(p*g*(g-1)),n-g)*sqrt(SIGMA_hat(3,3)*(1/n_time1+1/n_time2))
bonf_CI_mu_diff_13_X3 = mu_hat_time1(3)-mu_hat_time3(3) + [-1 1]*tinv(1-alpha/(p*g*(g-1)),n-g)*sqrt(SIGMA_hat(3,3)*(1/n_time1+1/n_time3))
bonf_CI_mu_diff_23_X3 = mu_hat_time2(3)-mu_hat_time3(3) + [-1 1]*tinv(1-alpha/(p*g*(g-1)),n-g)*sqrt(SIGMA_hat(3,3)*(1/n_time2+1/n_time3))
disp('-------------')
disp('Variable X4')
disp('-------------')
bonf_CI_mu_diff_12_X4 = mu_hat_time1(4)-mu_hat_time2(4) + [-1 1]*tinv(1-alpha/(p*g*(g-1)),n-g)*sqrt(SIGMA_hat(4,4)*(1/n_time1+1/n_time2))
bonf_CI_mu_diff_13_X4 = mu_hat_time1(4)-mu_hat_time3(4) + [-1 1]*tinv(1-alpha/(p*g*(g-1)),n-g)*sqrt(SIGMA_hat(4,4)*(1/n_time1+1/n_time3))
bonf_CI_mu_diff_23_X4 = mu_hat_time2(4)-mu_hat_time3(4) + [-1 1]*tinv(1-alpha/(p*g*(g-1)),n-g)*sqrt(SIGMA_hat(4,4)*(1/n_time2+1/n_time3))
disp('--------------------------------------------------------------------------')
%-------------------------------------------------------------------------------------------------------