clc;clear;close all
format compact
%--------------------------------------------------------------------------
mu = [0.5 1];
SIGMA = (1/12)*[1 0; 0 4];
p = 2;
n = 20;
runs = 10000;
X = zeros(runs,p);
for j = 1:runs
    Y = [rand(n,1) 2*rand(n,1)];
    X(j,:) = sqrt(n)*(mean(Y) - mu);
end
mu_X_hat = mean(X);
SIGMA_X_hat = cov(X);
%--------------------------------------------------------------------------
figure
subplot(5,2,1:2)
axis off
text(0,0.9,'Illustration of Central Limit Theorem for multivariate distributions (p = 2)','Fontsize',20)
text(0,0.3,'Y_i iid ~ Uniform2D([0,1],[0,2],  i = 1:n,  and Y_i_1,Y_i_2 independent    =>    {\mu}_Y = [1/2 1] = [0.5 1],  {\Sigma}_Y = (1/12)*[1 0; 0 4] = [0.0833 0; 0 0.333]','Fontsize',16) 
text(0,-0.05,'CLT:   Y  -->  N_2({\mu}_Y,{\Sigma}_Y/n)  for  n --> {\infty}    or equivalently    X = {\surd}n*(Y - {\mu}_Y)  -->  N_2(0,{\Sigma}_Y)  for  n --> {\infty}','Fontsize',16)
text(0.474,0.18,'_','Fontsize',12)
text(0.055,0.22,'_','Fontsize',16)
text(0.496,0.22,'_','Fontsize',16)

subplot(5,2,3)
axis off
x1_aver = mean(X(:,1));
x2_aver = mean(X(:,2));
S = cov(X);
text(0,0.75,'Descriptive statistics:','Fontsize',18)
text(0,1.1,['n = ' num2str(n)],'Fontsize',16)
text(0.45,0.75,[num2str(runs) ' runs  giving  X_j, j=1..' num2str(runs)],'Fontsize',16)
text(0.45,0.5,['{\mu}_X_,_h_a_t = x = ' mat2str([x1_aver x2_aver],2)],'Fontsize',16)
text(0.585,0.71,'_','Fontsize',16),warning('off')
text(0.45,0.25,['{\Sigma}_X_,_h_a_t = S_X = ' mat2str(S,3)],'Fontsize',16)

subplot(5,2,5:2:9)
plotmatrix(X,'r*')
xlabel('x_1                                          x_2','Fontsize',16)
ylabel('x_2                           x_1','Fontsize',16)
title(['Scattermatrix of X-observations for ' num2str(runs) ' runs'],'Fontsize',16)

subplot(5,2,4)
axis off
text(0,0.75,['Sample Mahalanobis distances, j=1..' num2str(runs) ':'],'Fontsize',18)
text(0,0.45,'d_j^2 = (x_j - x)^T S_X^-^1 (x_j - x) ~ {\chi}_2^2   (approximately)','Fontsize',16)
text(0.158,0.66,'_','Fontsize',16),warning('off')
text(0.358,0.66,'_','Fontsize',16),warning('off')

subplot(5,2,6:2:10)
mu_hat = [x1_aver x2_aver];
d_j_sqr = zeros(1,n);
for j = 1:runs
    d_j_sqr(j) = (X(j,:)-mu_hat)*inv(S)*(X(j,:)-mu_hat)';
end
df = 2;
z_j = chi2rnd(df,1,runs);
qqplot(d_j_sqr,z_j),grid,xlabel('quantiles for d_j^2','Fontsize',16),ylabel('quantiles for {\chi}_2^2 distribution','Fontsize',16),...
    title('qq-plot for d_j^2 versus {\chi}_2^2 distribution','Fontsize',16)
