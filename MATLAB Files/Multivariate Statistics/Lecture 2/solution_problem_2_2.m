clc;clear;close all
format compact
%--------------------------------------------------------------------------
disp('-----------------------------------------------------')
mu = [1; 2; 0]
SIGMA = [4 1 2; 1 9 -3; 2 -3 5]
[P,LAMBDA] = eig(SIGMA)
disp('-----------------------------------------------------')
P_LAMBDA_PT = P*LAMBDA*P'
disp('-----------------------------------------------------')
sqrtSIGMA = P*diag(sqrt(diag(LAMBDA)))*P'
command_sqrtmSIGMA = sqrtm(SIGMA)
disp('-----------------------------------------------------')
p = 3
n = 10000
disp('Yi ~ N(0,I3)')
Y = randn(n,p);
disp('Xi ~ mu + sqrtSIGMA*Yi')
%X = ones(n,1)*mu' + Y*sqrtSIGMA;
X = mu' + Y*sqrtSIGMA;
disp('-----------------------------------------------------')
mu_hat = mean(X)
SIGMA_hat = cov(X)
disp('-----------------------------------------------------')
%--------------------------------------------------------------------------
figure
subplot(5,2,1:2)
axis off
text(0,0.8,'Generation of general multivariate normal distribution (p = 3) from (univariate) standard normal distribution','Fontsize',20)
text(0,0.4,'Y_i ~ N_3(0,I_3)  =>  X_i = {\mu} + {\Sigma}^1^/^2*Y_i ~ N_3({\mu},{\Sigma}),  i=1..n','Fontsize',20)
text(0,0,['Moments:                        {\mu} = ' mat2str(mu') ',  {\Sigma} = ' mat2str(SIGMA)],'Fontsize',16)

subplot(5,2,3)
axis off
x1_aver = mean(X(:,1));
x2_aver = mean(X(:,2));
x3_aver = mean(X(:,3));
S = cov(X);
text(0,1,'Descriptive statistics:','Fontsize',18)
text(0.45,1,['n = ' num2str(n)],'Fontsize',16)
text(0.45,0.75,['{\mu}_h_a_t = x = ' mat2str([x1_aver x2_aver x3_aver],2)],'Fontsize',16)
text(0.555,0.96,'_','Fontsize',16),warning('off')
text(0.45,0.5,['{\Sigma}_h_a_t = S = ' mat2str(S,2)],'Fontsize',16)

subplot(5,2,5:2:9)
plotmatrix(X,'r*')
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
    d_i_sqr(i) = (X(i,:)-mu_hat)*inv(S)*(X(i,:)-mu_hat)';
end
df = 3;
z_i = chi2rnd(df,1,n);
qqplot(d_i_sqr,z_i),grid,xlabel('quantiles for d_i^2','Fontsize',16),ylabel('quantiles for {\chi}_3^2 distribution','Fontsize',16),...
    title('qq-plot for d_i^2 versus {\chi}_3^2 distribution','Fontsize',16)
%--------------------------------------------------------------------------
