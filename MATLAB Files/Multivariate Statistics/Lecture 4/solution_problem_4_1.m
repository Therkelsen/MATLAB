clc;clear;close all
format compact
%--------------------------------------------------------------------------
% load data
%--------------------------------------------------------------------------
load("Lecture 4/dataset_problem_4_1.mat")
disp('-------------------------------------------------------------------')
disp('Bartlett test,  H0: SIGMA_female = SIGMA_male')
disp('-------------------------------------------------------------------')
%--------------------------------------------------------------------------
% descriptive statistics
%--------------------------------------------------------------------------
p = 3;
n_female = 24;
n_male = 24;
Sigma_hat_female = cov(X_female)
Sigma_hat_male = cov(X_male)
Sigma_hat_pooled = ((n_female-1)*Sigma_hat_female + (n_male-1)*Sigma_hat_male)/(n_female+n_male-2)
%--------------------------------------------------------------------------
%  Bartlett test (Box's M-test) for equal covariance matrices
%--------------------------------------------------------------------------
T = (n_female+n_male-2)*log(det(Sigma_hat_pooled)) - (n_female-1)*log(det(Sigma_hat_female)) - (n_male-1)*log(det(Sigma_hat_male))
correction_factor = 1 - ((2*p^2+3*p-1)/(6*(p+1)))*(1/(n_female-1)+1/(n_male-1)-1/(n_female+n_male-2))
test_statistic = correction_factor*T
alpha = 0.05
df = p*(p+1)/2
critical_value = chi2inv(1-alpha,df)
reject_H0 =  test_statistic > critical_value
p_value = 1-chi2cdf(test_statistic,df)
disp('-------------------------------------------------------------------')
%------------------------------------------------------------------------------------------
% Model check for FEMALE data
%------------------------------------------------------------------------
figure
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for FEMALE observations possibly being from multivariate normal distribution (p = 3)','Fontsize',20)
subplot(5,2,3)
axis off
n = n_female;
x1_aver = mean(X_female(:,1));
x2_aver = mean(X_female(:,2));
x3_aver = mean(X_female(:,3));
S = cov(X_female);
text(0,1,'Descriptive statistics:','Fontsize',18)
text(0.45,1,['n = ' num2str(n)],'Fontsize',16)
text(0.45,0.75,['{\mu}_h_a_t = x = ' mat2str([x1_aver x2_aver x3_aver],2)],'Fontsize',16)
text(0.555,0.96,'_','Fontsize',16),warning('off')
text(0.45,0.5,['{\Sigma}_h_a_t = S = ' mat2str(S,2)],'Fontsize',16)
subplot(5,2,5:2:9)
plotmatrix(X_female,'ro')
xlabel('x_1                           x_2                          x_3','Fontsize',16)
ylabel('x_3                 x_2                 x_1','Fontsize',16)
title(['Scattermatrix of observations'],'Fontsize',16)
subplot(5,2,4)
axis off
text(0,1,'Sample Mahalanobis distances, i=1..n :','Fontsize',18)
text(0,0.7,'d_i^2 = (x_i - x)^T S^-^1 (x_i - x) ~ {\chi}_3^2   (approximately)','Fontsize',16)
text(0.158,0.91,'_','Fontsize',16),warning('off')
text(0.358,0.91,'_','Fontsize',16),warning('off')
subplot(5,2,6:2:10)
mu_hat = [x1_aver x2_aver x3_aver];
d_i_sqr = zeros(1,n);
for i = 1:n
    d_i_sqr(i) = (X_female(i,:)-mu_hat)*inv(S)*(X_female(i,:)-mu_hat)';
end
df = 3;
z_i = chi2rnd(df,1,n);
% z2 = chi2rnd(20,1,n);
qqplot(d_i_sqr,z_i),grid,xlabel('quantiles for d_i^2','Fontsize',16),ylabel('quantiles for {\chi}_3^2 distribution','Fontsize',16),...
    title('qq-plot for d_i^2 versus {\chi}_3^2 distribution','Fontsize',16)
%------------------------------------------------------------------------------------------
% Model check for MALE data
%------------------------------------------------------------------------
figure
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for MALE observations possibly being from multivariate normal distribution (p = 3)','Fontsize',20)
subplot(5,2,3)
axis off
n = n_male;
x1_aver = mean(X_male(:,1));
x2_aver = mean(X_male(:,2));
x3_aver = mean(X_male(:,3));
S = cov(X_male);
text(0,1,'Descriptive statistics:','Fontsize',18)
text(0.45,1,['n = ' num2str(n)],'Fontsize',16)
text(0.45,0.75,['{\mu}_h_a_t = x = ' mat2str([x1_aver x2_aver x3_aver],2)],'Fontsize',16)
text(0.555,0.96,'_','Fontsize',16),warning('off')
text(0.45,0.5,['{\Sigma}_h_a_t = S = ' mat2str(S,2)],'Fontsize',16)
subplot(5,2,5:2:9)
plotmatrix(X_male,'bo')
xlabel('x_1                           x_2                          x_3','Fontsize',16)
ylabel('x_3                 x_2                 x_1','Fontsize',16)
title(['Scattermatrix of observations'],'Fontsize',16)
subplot(5,2,4)
axis off
text(0,1,'Sample Mahalanobis distances, i=1..n :','Fontsize',18)
text(0,0.7,'d_i^2 = (x_i - x)^T S^-^1 (x_i - x) ~ {\chi}_3^2   (approximately)','Fontsize',16)
text(0.158,0.91,'_','Fontsize',16),warning('off')
text(0.358,0.91,'_','Fontsize',16),warning('off')
subplot(5,2,6:2:10)
mu_hat = [x1_aver x2_aver x3_aver];
d_i_sqr = zeros(1,n);
for i = 1:n
    d_i_sqr(i) = (X_male(i,:)-mu_hat)*inv(S)*(X_male(i,:)-mu_hat)';
end
df = 3;
z_i = chi2rnd(df,1,n);
% z2 = chi2rnd(20,1,n);
qqplot(d_i_sqr,z_i),grid,xlabel('quantiles for d_i^2','Fontsize',16),ylabel('quantiles for {\chi}_3^2 distribution','Fontsize',16),...
    title('qq-plot for d_i^2 versus {\chi}_3^2 distribution','Fontsize',16)
%------------------------------------------------------------------------