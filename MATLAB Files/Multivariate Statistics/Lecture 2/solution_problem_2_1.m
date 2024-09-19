clc;clear;close all
format compact

load('Lecture 2\dataset_problem_2_1.mat');
n = 130

figure(1)
subplot(5,2,1:2)
axis off
text(0,0.6,'Model check for observations possibly being from multivariate normal distribution (p = 5)','Fontsize',20)

subplot(5,2,3)
axis off
mu_hat = mean(X)
S = cov(X)
R = corrcoef(X)
text(0,1,'Descriptive statistics:','Fontsize',18)
text(0,0.7,['n = ' num2str(n)],'Fontsize',16)
text(0,0.4,'see command window for  {\mu}_h_a_t = x,  {\Sigma}_h_a_t = S  and  {\rho}_h_a_t = R','Fontsize',16)
text(0.59,0.63,'_','Fontsize',16)

subplot(5,2,5:2:9)
plotmatrix(X,'r*')
xlabel('x_1              x_2              x_3              x_4              x_5','Fontsize',16)
ylabel('x_5       x_4         x_3         x_2         x_1','Fontsize',16)
title(['Scattermatrix of observations'],'Fontsize',16)

subplot(5,2,4)
axis off
text(0,1,'Sample Mahalanobis distances, i=1..n :','Fontsize',18)
text(0,0.7,'d_i^2 = (x_i - x)^T S^-^1 (x_i - x) ~ {\chi}_5^2   (approximately)','Fontsize',16)
text(0.158,0.91,'_','Fontsize',16),warning('off')
text(0.358,0.91,'_','Fontsize',16),warning('off')

subplot(5,2,6:2:10)
d_i_sqr = zeros(1,n);
for i = 1:n
    d_i_sqr(i) = (X(i,:)-mu_hat)*inv(S)*(X(i,:)-mu_hat)';
end
df = 5;
z_i = chi2rnd(df,1,10000);
qqplot(d_i_sqr,z_i),grid,xlabel('quantiles for d_i^2','Fontsize',16),ylabel('quantiles for {\chi}_5^2 distribution','Fontsize',16),...
    title('qq-plot for d_i^2 versus {\chi}_5^2 distribution','Fontsize',16)

%-----------------------------------------------------------------------------------------------------------------
