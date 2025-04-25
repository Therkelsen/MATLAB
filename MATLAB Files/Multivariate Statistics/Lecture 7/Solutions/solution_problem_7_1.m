clc;clear;close all
format compact
%--------------------------------------------------------------------------
% load data
%--------------------------------------------------------------------------
% manual calculation
%--------------------------------------------------------------------------
load('C:\Users\cv\OneDrive - Syddansk Universitet\M-drive moved by SDU-IT\Statistik\MULTISTAT\Lektion 07  MMLR\Course material\dataset_problem_7_1');
disp('-------------------------------------------------------------------')
disp('Multivariate Multiple regression model:  Y = Zb + e')
disp('                                         Y = (n x m)')
disp('                                         Z = (n x (r+1))')
disp('                                         b = ((r+1) x m)')
disp('                                         e = (n x m)')
disp('-------------------------------------------------------------------')
n = size(Y,1)
m = size(Y,2)
r = size(Z,2)
Z = [ones(n,1) Z];
disp('-------------------------------------------------------------------')
disp('Unbiased Least Squares estimate of regression coefficients')
disp('-------------------------------------------------------------------')
b_hat = inv(Z'*Z)*Z'*Y
disp('-------------------------------------------------------------------')
disp('Unbiased estimate of SIGMA_e = Cov[e]        (m x m)')
disp('-------------------------------------------------------------------')
Y_hat = Z*b_hat;
e_hat = Y - Y_hat;
SSE = e_hat'*e_hat;
SIGMA_e_hat = SSE/(n-(r+1))
test = cov(e_hat)*(n-1)/(n-r-1)
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp('Uncorrected SS decomposition: SST = SSR + SSE     (m x m)')
disp('-------------------------------------------------------------------')
SST_uncorrected = Y'*Y
SSR_uncorrected = Y_hat'*Y_hat
SSE 
disp('-------------------------------------------------------------------')
disp('Mean corrected SS decomposition: SST = SSR + SSE   (m x m)')
disp('-------------------------------------------------------------------')
SST_mean_corrected = (Y-ones(n,1)*mean(Y))'*(Y-ones(n,1)*mean(Y))
SSR_mean_corrected = (Y_hat-ones(n,1)*mean(Y))'*(Y_hat-ones(n,1)*mean(Y))
SSE
disp('-------------------------------------------------------------------')
disp('Coefficient of determination, R_sqr, and R_adj_sqr')
disp('-------------------------------------------------------------------')
R_sqr = 1 - det(SSE)/det(SST_mean_corrected)
R_adj_sqr = 1 - (det(SSE)/(n-(r+1)*m))/(det(SST_mean_corrected)/(n-m))
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp('Inference for regression coeffs for Y1,  b_x1 = [b_01 ... b_r1]')
disp('-------------------------------------------------------------------')
disp('Estimate of b_x1:')
disp('-------------------------------------------------------------------')
b_x1_hat = b_hat(:,1)
disp('-------------------------------------------------------------------')
disp('Estimated covariance matrix for b_x1_hat   ((r+1) x (r+1))')
disp('-------------------------------------------------------------------')
SIGMA_hat_b_x1_hat = SIGMA_e_hat(1,1)*inv(Z'*Z)
sigma_hat_b_01_hat = sqrt(SIGMA_hat_b_x1_hat(1,1))
sigma_hat_b_11_hat = sqrt(SIGMA_hat_b_x1_hat(2,2))
sigma_hat_b_21_hat = sqrt(SIGMA_hat_b_x1_hat(3,3))
sigma_hat_b_31_hat = sqrt(SIGMA_hat_b_x1_hat(4,4))
sigma_hat_b_41_hat = sqrt(SIGMA_hat_b_x1_hat(5,5))
disp('-------------------------------------------------------------------')
disp('95% simultaneous confidence intervals for b_x1')
disp('-------------------------------------------------------------------')
alpha = 0.05;
sim_CI_b_01 = b_x1_hat(1) + [-1 1]*sqrt((r+1)*finv(1-alpha,r+1,n-(r+1)))*sigma_hat_b_01_hat
sim_CI_b_11 = b_x1_hat(2) + [-1 1]*sqrt((r+1)*finv(1-alpha,r+1,n-(r+1)))*sigma_hat_b_11_hat
sim_CI_b_21 = b_x1_hat(3) + [-1 1]*sqrt((r+1)*finv(1-alpha,r+1,n-(r+1)))*sigma_hat_b_21_hat
sim_CI_b_31 = b_x1_hat(4) + [-1 1]*sqrt((r+1)*finv(1-alpha,r+1,n-(r+1)))*sigma_hat_b_31_hat
sim_CI_b_41 = b_x1_hat(5) + [-1 1]*sqrt((r+1)*finv(1-alpha,r+1,n-(r+1)))*sigma_hat_b_41_hat
disp('-------------------------------------------------------------------')
disp('95% Bonferroni confidence intervals for b_x1')
disp('-------------------------------------------------------------------')
bonf_CI_b_01 = b_x1_hat(1) + [-1 1]*tinv(1-alpha/(2*(r+1)),n-(r+1))*sigma_hat_b_01_hat
bonf_CI_b_11 = b_x1_hat(2) + [-1 1]*tinv(1-alpha/(2*(r+1)),n-(r+1))*sigma_hat_b_11_hat
bonf_CI_b_21 = b_x1_hat(3) + [-1 1]*tinv(1-alpha/(2*(r+1)),n-(r+1))*sigma_hat_b_21_hat
bonf_CI_b_31 = b_x1_hat(4) + [-1 1]*tinv(1-alpha/(2*(r+1)),n-(r+1))*sigma_hat_b_31_hat
bonf_CI_b_41 = b_x1_hat(5) + [-1 1]*tinv(1-alpha/(2*(r+1)),n-(r+1))*sigma_hat_b_41_hat
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp('Inference for regression coeffs for Y2,  b_x2 = [b_02 ... b_r2])')
disp('Inference for regression coeffs for Y3,  b_x3 = [b_03 ... b_r3])')
disp('Inference for regression coeffs for Y4,  b_x4 = [b_04 ... b_r4])')
disp('(analogue to inference for regress.coeffs for Y1 and NOT performed here)')
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp('Extra Sum of Squares based test for significant regression model') 
disp('(significant dependency of Y = [Y1 ... Ym] on [z1 z2 ... zr])')
disp('H0: all regress.coeffs = 0 except intercept coeffs b_01,...,b_0m')
disp('-------------------------------------------------------------------')
q = 0
SSE_full_model = (Y - Z*b_hat)'*(Y - Z*b_hat)
Z_only_intercept = Z(:,1:(q+1));
b_hat_only_intercept = inv(Z_only_intercept'*Z_only_intercept)*Z_only_intercept'*Y
SSE_only_intercept_model = (Y - Z_only_intercept*b_hat_only_intercept)'*(Y - Z_only_intercept*b_hat_only_intercept)
LAMBDA_LR = (det(SSE_full_model/n)/det(SSE_only_intercept_model/n))^(n/2)
test_statistic = -(n-r-1-(m-r+q+1)/2)*log(LAMBDA_LR^(2/n))
critical_value = chi2inv(1-alpha,m*(r-q)) 
rejection_of_H0 = test_statistic >= critical_value
p_value = 1-chi2cdf(test_statistic,m*(r-q))
disp('-------------------------------------------------------------------')
disp('Extra Sum of Squares based test for significant predictor variable z1') 
disp('(significant dependency of Y = [Y1 ... Ym] on z1)')
disp('H0: b_11 = ... = b_1m = 0')
disp('-------------------------------------------------------------------')
q = 3
SSE_full_model = (Y - Z*b_hat)'*(Y - Z*b_hat)
Z_model_without_z1 = Z(:,[1 3:(r+1)]);
b_hat_model_without_z1 = inv(Z_model_without_z1'*Z_model_without_z1)*Z_model_without_z1'*Y
SSE_model_without_z1 = (Y - Z_model_without_z1*b_hat_model_without_z1)'*(Y - Z_model_without_z1*b_hat_model_without_z1)
LAMBDA_LR = (det(SSE_full_model/n)/det(SSE_model_without_z1/n))^(n/2)
test_statistic = -(n-r-1-(m-r+q+1)/2)*log(LAMBDA_LR^(2/n))
critical_value = chi2inv(1-alpha,m*(r-q)) 
rejection_of_H0 = test_statistic >= critical_value
p_value = 1-chi2cdf(test_statistic,m*(r-q))
disp('-------------------------------------------------------------------')
disp('Extra Sum of Squares based test for significant predictor variable z2') 
disp('(significant dependency of Y = [Y1 ... Ym] on z2)')
disp('H0: b_21 = ... = b_2m = 0')
disp('-------------------------------------------------------------------')
q = 3
SSE_full_model = (Y - Z*b_hat)'*(Y - Z*b_hat)
Z_model_without_z2 = Z(:,[1:2 4:(r+1)]);
b_hat_model_without_z2 = inv(Z_model_without_z2'*Z_model_without_z2)*Z_model_without_z2'*Y
SSE_model_without_z2 = (Y - Z_model_without_z2*b_hat_model_without_z2)'*(Y - Z_model_without_z2*b_hat_model_without_z2)
LAMBDA_LR = (det(SSE_full_model/n)/det(SSE_model_without_z2/n))^(n/2)
test_statistic = -(n-r-1-(m-r+q+1)/2)*log(LAMBDA_LR^(2/n))
critical_value = chi2inv(1-alpha,m*(r-q)) 
rejection_of_H0 = test_statistic >= critical_value
p_value = 1-chi2cdf(test_statistic,m*(r-q))
disp('-------------------------------------------------------------------')
disp('Extra Sum of Squares based test for significant predictor variable z3') 
disp('(significant dependency of Y = [Y1 ... Ym] on z3)')
disp('H0: b_31 = ... = b_3m = 0')
disp('-------------------------------------------------------------------')
q = 3
SSE_full_model = (Y - Z*b_hat)'*(Y - Z*b_hat)
Z_model_without_z3 = Z(:,[1:3 5:(r+1)]);
b_hat_model_without_z3 = inv(Z_model_without_z3'*Z_model_without_z3)*Z_model_without_z3'*Y
SSE_model_without_z3 = (Y - Z_model_without_z3*b_hat_model_without_z3)'*(Y - Z_model_without_z3*b_hat_model_without_z3)
LAMBDA_LR = (det(SSE_full_model/n)/det(SSE_model_without_z3/n))^(n/2)
test_statistic = -(n-r-1-(m-r+q+1)/2)*log(LAMBDA_LR^(2/n))
critical_value = chi2inv(1-alpha,m*(r-q)) 
rejection_of_H0 = test_statistic >= critical_value
p_value = 1-chi2cdf(test_statistic,m*(r-q))
disp('-------------------------------------------------------------------')
disp('Extra Sum of Squares based test for significant predictor variable z4') 
disp('(significant dependency of Y = [Y1 ... Ym] on z4)')
disp('H0: b_41 = ... = b_4m = 0')
disp('-------------------------------------------------------------------')
q = 3
SSE_full_model = (Y - Z*b_hat)'*(Y - Z*b_hat)
Z_model_without_z4 = Z(:,1:r);
b_hat_model_without_z4 = inv(Z_model_without_z4'*Z_model_without_z4)*Z_model_without_z4'*Y
SSE_model_without_z4 = (Y - Z_model_without_z4*b_hat_model_without_z4)'*(Y - Z_model_without_z4*b_hat_model_without_z4)
LAMBDA_LR = (det(SSE_full_model/n)/det(SSE_model_without_z4/n))^(n/2)
test_statistic = -(n-r-1-(m-r+q+1)/2)*log(LAMBDA_LR^(2/n))
critical_value = chi2inv(1-alpha,m*(r-q)) 
rejection_of_H0 = test_statistic >= critical_value
p_value = 1-chi2cdf(test_statistic,m*(r-q))
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp('95% simultaneous prediction intervals for E[Y(z0)]')
disp('(mean of new responses at z=z0)')
disp('-------------------------------------------------------------------')
z0 = [1 0.330 45.500 20.375 1.010]
Y_z0_hat = z0*b_hat
sim_pred_CI_mean_Y1_z0 = Y_z0_hat(1) + [-1 0 1]*sqrt(z0*inv(Z'*Z)*z0'*(n/(n-r-1))*SIGMA_e_hat(1,1))*sqrt((m*(n-r-1)/(n-r-m))*finv(1-alpha,m,n-r-m))
sim_pred_CI_mean_Y2_z0 = Y_z0_hat(2) + [-1 0 1]*sqrt(z0*inv(Z'*Z)*z0'*(n/(n-r-1))*SIGMA_e_hat(2,2))*sqrt((m*(n-r-1)/(n-r-m))*finv(1-alpha,m,n-r-m))
sim_pred_CI_mean_Y3_z0 = Y_z0_hat(3) + [-1 0 1]*sqrt(z0*inv(Z'*Z)*z0'*(n/(n-r-1))*SIGMA_e_hat(3,3))*sqrt((m*(n-r-1)/(n-r-m))*finv(1-alpha,m,n-r-m))
sim_pred_CI_mean_Y4_z0 = Y_z0_hat(4) + [-1 0 1]*sqrt(z0*inv(Z'*Z)*z0'*(n/(n-r-1))*SIGMA_e_hat(4,4))*sqrt((m*(n-r-1)/(n-r-m))*finv(1-alpha,m,n-r-m))
disp('-------------------------------------------------------------------')
disp('95% prediction interval for Y(z0)')
disp('(single new response at z=z0)')
disp('-------------------------------------------------------------------')
z0 = [1 0.330 45.500 20.375 1.010]
Y_z0_hat = z0*b_hat
sim_pred_CI_Y1_z0 = Y_z0_hat(1) + [-1 0 1]*sqrt((1+z0*inv(Z'*Z)*z0')*(n/(n-r-1))*SIGMA_e_hat(1,1))*sqrt((m*(n-r-1)/(n-r-m))*finv(1-alpha,m,n-r-m))
sim_pred_CI_Y2_z0 = Y_z0_hat(2) + [-1 0 1]*sqrt((1+z0*inv(Z'*Z)*z0')*(n/(n-r-1))*SIGMA_e_hat(2,2))*sqrt((m*(n-r-1)/(n-r-m))*finv(1-alpha,m,n-r-m))
sim_pred_CI_Y3_z0 = Y_z0_hat(3) + [-1 0 1]*sqrt((1+z0*inv(Z'*Z)*z0')*(n/(n-r-1))*SIGMA_e_hat(3,3))*sqrt((m*(n-r-1)/(n-r-m))*finv(1-alpha,m,n-r-m))
sim_pred_CI_Y4_z0 = Y_z0_hat(4) + [-1 0 1]*sqrt((1+z0*inv(Z'*Z)*z0')*(n/(n-r-1))*SIGMA_e_hat(4,4))*sqrt((m*(n-r-1)/(n-r-m))*finv(1-alpha,m,n-r-m))
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp('Model check: univariate analysis of residuals (estimated model errors)')
disp('SEE FIGURE 1-4 for marginal analysis of e_j1,...,e_j4')
disp('-------------------------------------------------------------------')
e1_hat = Y(:,1) - Y_hat(:,1);
figure(1)
subplot(7,3,1:3)
axis off
text(0.04,0.4,'Analysis of residuals (estimated model errors),  e_j_1_,_h_a_t = Y_j_1 - Y_j_1_,_h_a_t = Y_j_1 - Z*b_h_a_t,  j = 1,2,...,n','Fontsize',20)
subplot(7,3,[7 10])
j = 1:n;
plot(j,e1_hat,'r+'),line([1 n],[0 0],'color','k'),xlim([1 n]),title('e_j_1_,_h_a_t = f(j)','Fontsize',14)
subplot(7,3,[8 11])
plot(Y(:,1),e1_hat,'k+'),line([min(Y(:,1)) max(Y(:,1))],[0 0],'color','k'),xlim([min(Y(:,1)) max(Y(:,1))]),title('e_j_1_,_h_a_t = f(Y_j_1)','Fontsize',14)
subplot(7,3,[9 12])
hist(e1_hat),title('histogram of e_j_1_,_h_a_t','Fontsize',14)
subplot(7,3,[18 21])
normplot(e1_hat),title('normal probability plot (QQ-plot)','Fontsize',14)
%--------------------------------------------------------------------------
e2_hat = Y(:,2) - Y_hat(:,2);
figure(2)
subplot(7,3,1:3)
axis off
text(0.04,0.4,'Analysis of residuals (estimated model errors),  e_j_2_,_h_a_t = Y_j_2 - Y_j_2_,_h_a_t = Y_j_2 - Z*b_h_a_t,  j = 1,2,...,n','Fontsize',20)
subplot(7,3,[7 10])
j = 1:n;
plot(j,e2_hat,'r+'),line([1 n],[0 0],'color','k'),xlim([1 n]),title('e_j_2_,_h_a_t = f(j)','Fontsize',14)
subplot(7,3,[8 11])
plot(Y(:,2),e2_hat,'k+'),line([min(Y(:,2)) max(Y(:,2))],[0 0],'color','k'),xlim([min(Y(:,2)) max(Y(:,2))]),title('e_j_2_,_h_a_t = f(Y_j_2)','Fontsize',14)
subplot(7,3,[9 12])
hist(e2_hat),title('histogram of e_j_2_,_h_a_t','Fontsize',14)
subplot(7,3,[18 21])
normplot(e2_hat),title('normal probability plot (QQ-plot)','Fontsize',14)
%--------------------------------------------------------------------------
e3_hat = Y(:,3) - Y_hat(:,3);
figure(3)
subplot(7,3,1:3)
axis off
text(0.04,0.4,'Analysis of residuals (estimated model errors),  e_j_3_,_h_a_t = Y_j_3 - Y_j_3_,_h_a_t = Y_j_3 - Z*b_h_a_t,  j = 1,2,...,n','Fontsize',20)
subplot(7,3,[7 10])
j = 1:n;
plot(j,e3_hat,'r+'),line([1 n],[0 0],'color','k'),xlim([1 n]),title('e_j_3_,_h_a_t = f(j)','Fontsize',14)
subplot(7,3,[8 11])
plot(Y(:,3),e3_hat,'k+'),line([min(Y(:,3)) max(Y(:,3))],[0 0],'color','k'),xlim([min(Y(:,3)) max(Y(:,3))]),title('e_j_3_,_h_a_t = f(Y_j_3)','Fontsize',14)
subplot(7,3,[9 12])
hist(e3_hat),title('histogram of e_j_3_,_h_a_t','Fontsize',14)
subplot(7,3,[18 21])
normplot(e3_hat),title('normal probability plot (QQ-plot)','Fontsize',14)
%--------------------------------------------------------------------------
e4_hat = Y(:,4) - Y_hat(:,4);
figure(4)
subplot(7,3,1:3)
axis off
text(0.04,0.4,'Analysis of residuals (estimated model errors),  e_j_4_,_h_a_t = Y_j_4 - Y_j_4_,_h_a_t = Y_j_4 - Z*b_h_a_t,  j = 1,2,...,n','Fontsize',20)
subplot(7,3,[7 10])
j = 1:n;
plot(j,e4_hat,'r+'),line([1 n],[0 0],'color','k'),xlim([1 n]),title('e_j_4_,_h_a_t = f(j)','Fontsize',14)
subplot(7,3,[8 11])
plot(Y(:,4),e4_hat,'k+'),line([min(Y(:,4)) max(Y(:,4))],[0 0],'color','k'),xlim([min(Y(:,4)) max(Y(:,4))]),title('e_j_4_,_h_a_t = f(Y_j_4)','Fontsize',14)
subplot(7,3,[9 12])
hist(e4_hat),title('histogram of e_j_4_,_h_a_t','Fontsize',14)
subplot(7,3,[18 21])
normplot(e4_hat),title('normal probability plot (QQ-plot)','Fontsize',14)
disp('-------------------------------------------------------------------')
disp('Model check: multivariate analysis of residuals (estimated model errors)')
disp('SEE FIGURE 5 for simultaneous analysis of [e_j1,...,e_j4]')
disp('-------------------------------------------------------------------')
e_hat = Y - Y_hat;
figure(5)
subplot(5,2,1:2)
axis off
text(0,0.3,'Model check for residuals possibly being from multivariate normal distribution (p = 4)','Fontsize',20)
subplot(5,2,5:2:9)
plotmatrix(e_hat,'r*')
title('Scattermatrix of multivariate residuals','Fontsize',16)
subplot(5,2,4)
axis off
text(0,0.5,'Sample Mahalanobis distances, j=1..n :','Fontsize',18)
subplot(5,2,6:2:10)
d_j_sqr = zeros(1,n);
for j = 1:n
    d_j_sqr(j) = (e_hat(j,:)-mean(e_hat))*inv(SIGMA_e_hat)*(e_hat(j,:)-mean(e_hat))';
end
df = 4;
z_j = chi2rnd(df,1,n);
qqplot(d_j_sqr,z_j),grid,xlabel('quantiles for d_j^2','Fontsize',16),ylabel('quantiles for {\chi}_4^2 distribution','Fontsize',16),...
    title('qq-plot for d_j^2 versus {\chi}_4^2 distribution','Fontsize',16)
%--------------------------------------------------------------------------
% calculations using built-in functions
%--------------------------------------------------------------------------
disp('-------------------------------------------------------------------')
disp('Using built-in function:  mvregress')
disp('-------------------------------------------------------------------')
Z_cells = cell(n,1);
for j = 1:n
	Z_cells{j} = [];
    for k = 1:(r+1)
        Z_cells{j} = [Z_cells{j} Z(j,k)*eye(m)];
    end
end
[beta,Sigma,E,CovB,logL] = mvregress(Z_cells,Y);
b_hat_mvregress = (reshape(beta,m,r+1))'
b_hat
disp('-------------------------------------------------------------------')
SIGMA_e_hat_mvregress = (n/(n-(r+1)))*cov(E)
SIGMA_e_hat
disp('-------------------------------------------------------------------')
% %-------------------------------------------------------------------------------------------------------