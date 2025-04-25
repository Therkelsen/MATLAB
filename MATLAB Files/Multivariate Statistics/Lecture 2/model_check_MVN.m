clc;clear;close all
format compact

%------------------------------------------------------------------------
% p = 1 with fitting model
%------------------------------------------------------------------------
mu = 10;
sigma = 2;
n = 1000;
x = mu + sigma*randn(1,n);

figure(1)
subplot(5,3,1:3)
axis off
text(0,0.7,'Model check for observations possibly being from univariate normal distribution (~ MVN with p = 1)','Fontsize',20)

subplot(5,3,4:3:13)
axis off
x_aver = mean(x);
Sxx = var(x);
text(0,0.9,'Descriptive statistics:','Fontsize',18)
text(0,0.8,['n = ' num2str(n)],'Fontsize',16)
text(0,0.7,['x = ' num2str(x_aver,2)],'Fontsize',16)
text(0,0.74,'_','Fontsize',16),warning('off')
text(0,0.6,['S_x_x = S^2 =  ' num2str(Sxx,3)],'Fontsize',16)
text(0,0.5,['s =  ' num2str(sqrt(Sxx),3)],'Fontsize',16)

subplot(5,3,5:3:14)
hist(x,25),title('histogram for observations','Fontsize',14),xlim([3 17]),xlabel('x-intervals (bins)','Fontsize',14)

subplot(5,3,6:3:15)
normplot(x),title('normal probability plot (qq-plot)','Fontsize',14),xlabel('data values, x_i','Fontsize',14),ylabel('quantiles for data (blue) and for normal distribution (red)','Fontsize',14)

%------------------------------------------------------------------------
% p = 1 with non-fitting model
%------------------------------------------------------------------------
df = 10;
n = 1000;
x = chi2rnd(df,1,n);

figure(2)
subplot(5,3,1:3)
axis off
text(0,0.7,'Model check for observations possibly being from univariate normal distribution (~ MVN with p = 1)','Fontsize',20)

subplot(5,3,4:3:7)
axis off
x_aver = mean(x);
Sxx = var(x);
text(0,0.9,'Descriptive statistics:','Fontsize',18)
text(0,0.75,['n = ' num2str(n)],'Fontsize',16)
text(0,0.6,['x = ' num2str(x_aver,2)],'Fontsize',16)
text(0,0.68,'_','Fontsize',16),warning('off')
text(0,0.45,['S_x_x = S_x^2 =  ' num2str(Sxx,3)],'Fontsize',16)
text(0,0.3,['s_x =  ' num2str(sqrt(Sxx),3)],'Fontsize',16)

subplot(5,3,5:3:8)
hist(x,25),title('histogram for observations','Fontsize',14),xlim([0 30])

subplot(5,3,6:3:9)
normplot(x),title('qq-plot for observations','Fontsize',14),xlabel(' '),ylabel('quantiles','Fontsize',14)

y = sqrt(x);

subplot(5,3,10:3:13)
axis off
y_aver = mean(y);
Syy = var(y);
text(0,1.15,'Transformation of observations:','Fontsize',20)
text(0,1,'y_i = sqrt(x_i),  i = 1..n','Fontsize',20)

text(0,0.7,'Descriptive statistics:','Fontsize',18)
text(0,0.55,['y = ' num2str(y_aver,2)],'Fontsize',16)
text(0,0.63,'_','Fontsize',16),warning('off')
text(0,0.4,['S_y_y = S_y^2 =  ' num2str(Syy,3)],'Fontsize',16)
text(0,0.25,['s_y =  ' num2str(sqrt(Syy),3)],'Fontsize',16)

subplot(5,3,11:3:14)
hist(y,25),xlim([0 7]),xlabel('histogram for transformed observations','Fontsize',14)

subplot(5,3,12:3:15)
normplot(y),title(' '),xlabel('qq-plot for transformed observations','Fontsize',14),ylabel('quantiles','Fontsize',14)

%------------------------------------------------------------------------
% p = 2 with fitting model
%------------------------------------------------------------------------
mu = [0 0];
sigma1 = 2;
sigma2 = 4;
rho12 = 0.8;
SIGMA = [sigma1^2 rho12*sigma1*sigma2; rho12*sigma1*sigma2 sigma2^2];
n = 1000;
X = mvnrnd(mu,SIGMA,n);

figure(3)
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for observations possibly being from bivariate normal distribution (~ MVN with p = 2)','Fontsize',20)

subplot(5,2,3)
axis off
x1_aver = mean(X(:,1));
x2_aver = mean(X(:,2));
S = cov(X);
text(0,1,'Descriptive statistics:','Fontsize',18)
text(0.5,1,['n = ' num2str(n)],'Fontsize',16)
text(0.5,0.75,['{\mu}_h_a_t = x = ' mat2str([x1_aver x2_aver],2)],'Fontsize',16)
text(0.605,0.96,'_','Fontsize',16),warning('off')
text(0.5,0.5,['{\Sigma}_h_a_t = S = ' mat2str(S,3)],'Fontsize',16)

subplot(5,2,5:2:9)
plotmatrix(X,'r*')
xlabel('x_1                                          x_2','Fontsize',16)
ylabel('x_2                           x_1','Fontsize',16)
title(['Scattermatrix of observations'],'Fontsize',16)

subplot(5,2,4)
axis off
text(0,1,'Sample Mahalanobis distances, i=1..n :','Fontsize',18)
text(0,0.7,'d_i^2 = (x_i - x)^T S^-^1 (x_i - x) ~ {\chi}_2^2   (approximately)','Fontsize',16)
text(0.158,0.91,'_','Fontsize',16),warning('off')
text(0.358,0.91,'_','Fontsize',16),warning('off')

subplot(5,2,6:2:10)
mu_hat = [x1_aver x2_aver];
d_i_sqr = zeros(1,n);
for i = 1:n
    d_i_sqr(i) = (X(i,:)-mu_hat)*inv(S)*(X(i,:)-mu_hat)';
end
df = 2;
z_i = chi2rnd(df,1,n);
% z2 = chi2rnd(20,1,n);
qqplot(d_i_sqr,z_i),grid,xlabel('quantiles for d_i^2','Fontsize',16),ylabel('quantiles for {\chi}_2^2 distribution','Fontsize',16),...
    title('qq-plot for d_i^2 versus {\chi}_2^2 distribution','Fontsize',16)

%------------------------------------------------------------------------
% p = 2 with non-fitting model
%------------------------------------------------------------------------
mu = [0 0];
sigma1 = 2;
sigma2 = 4;
rho12 = 0.8;
SIGMA = [sigma1^2 rho12*sigma1*sigma2; rho12*sigma1*sigma2 sigma2^2];
n = 1000;
X = mvnrnd(mu,SIGMA,n);

X(:,1) = X(:,1).^2;

figure(4)
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for observations possibly being from bivariate normal distribution (~ MVN with p = 2)','Fontsize',20)

subplot(5,2,3)
axis off
x1_aver = mean(X(:,1));
x2_aver = mean(X(:,2));
S = cov(X);
text(0,1,'Descriptive statistics:','Fontsize',18)
text(0.5,1,['n = ' num2str(n)],'Fontsize',16)
text(0.5,0.75,['{\mu}_h_a_t = x = ' mat2str([x1_aver x2_aver],2)],'Fontsize',16)
text(0.605,0.96,'_','Fontsize',16),warning('off')
text(0.5,0.5,['{\Sigma}_h_a_t = S = ' mat2str(S,3)],'Fontsize',16)

subplot(5,2,5:2:9)
plotmatrix(X,'r*')
xlabel('x_1                                          x_2','Fontsize',16)
ylabel('x_2                           x_1','Fontsize',16)
title(['Scattermatrix of observations'],'Fontsize',16)

subplot(5,2,4)
axis off
text(0,1,'Sample Mahalanobis distances, i=1..n :','Fontsize',18)
text(0,0.7,'d_i^2 = (x_i - x)^T S^-^1 (x_i - x) ~ {\chi}_2^2   (approximately)','Fontsize',16)
text(0.158,0.91,'_','Fontsize',16),warning('off')
text(0.358,0.91,'_','Fontsize',16),warning('off')

subplot(5,2,6:2:10)
mu_hat = [x1_aver x2_aver];
d_i_sqr = zeros(1,n);
for i = 1:n
    d_i_sqr(i) = (X(i,:)-mu_hat)*inv(S)*(X(i,:)-mu_hat)';
end
df = 2;
z_i = chi2rnd(df,1,n);
% z2 = chi2rnd(20,1,n);
qqplot(d_i_sqr,z_i),grid,xlabel('quantiles for d_i^2','Fontsize',16),ylabel('quantiles for {\chi}_2^2 distribution','Fontsize',16),...
    title('qq-plot for d_i^2 versus {\chi}_2^2 distribution','Fontsize',16)
 
%------------------------------------------------------------------------
% p = 3 with fitting model
%------------------------------------------------------------------------
mu = [0 0 0];
sigma1 = 2;
sigma2 = 4;
sigma3 = 2;
rho12 = 0.95;
rho13 = -0.7;
rho23 = -0.5;
SIGMA = [sigma1^2 rho12*sigma1*sigma2 rho13*sigma1*sigma3; rho12*sigma1*sigma2 sigma2^2 rho23*sigma2*sigma3; rho13*sigma1*sigma3 rho23*sigma2*sigma3 sigma3^2];
n = 1000;
X = mvnrnd(mu,SIGMA,n);

figure(5)
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for observations possibly being from multivariate normal distribution (p = 3)','Fontsize',20)

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
% z2 = chi2rnd(20,1,n);
qqplot(d_i_sqr,z_i),grid,xlabel('quantiles for d_i^2','Fontsize',16),ylabel('quantiles for {\chi}_3^2 distribution','Fontsize',16),...
    title('qq-plot for d_i^2 versus {\chi}_3^2 distribution','Fontsize',16)

%------------------------------------------------------------------------
% p = 3 with non-fitting model
%------------------------------------------------------------------------
mu = [0 0 0];
sigma1 = 2;
sigma2 = 4;
sigma3 = 2;
rho12 = 0.95;
rho13 = -0.7;
rho23 = -0.5;
SIGMA = [sigma1^2 rho12*sigma1*sigma2 rho13*sigma1*sigma3; rho12*sigma1*sigma2 sigma2^2 rho23*sigma2*sigma3; rho13*sigma1*sigma3 rho23*sigma2*sigma3 sigma3^2];
n = 1000;
X = mvnrnd(mu,SIGMA,n);

X(:,2) = X(:,2).^2;

figure(6)
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for observations possibly being from multivariate normal distribution (p = 3)','Fontsize',20)

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
test = mahal(X,X);
df = 3;
z_i = chi2rnd(df,1,n);
% z2 = chi2rnd(20,1,n);
qqplot(d_i_sqr,z_i),grid,xlabel('quantiles for d_i^2','Fontsize',16),ylabel('quantiles for {\chi}_3^2 distribution','Fontsize',16),...
    title('qq-plot for d_i^2 versus {\chi}_3^2 distribution','Fontsize',16)
%--------------------------------------------------------------------------