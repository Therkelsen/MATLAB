clc;clear;close all
format compact 

subplot(5,1,1)
axis off
text(0,1.2,'Hypothesis testing for one-sample multivariate normal population mean vector,  H_o: {\mu} = {\mu}_o','Fontsize',16)
text(0,0.8,'One-to-one, monotone equivalence between Hotelling test-statistic,                  T^2 = (x - {\mu}_o)^T (S/n)^-^1 (x - {\mu}_o),       0 <= T^2 < {\infty}','Fontsize',14)
text(0.572,0.99,'_','Fontsize',16),warning('off')
text(0.682,0.99,'_','Fontsize',16),warning('off')
text(0.275,0.55,'and Likelihood Ratio test-statistic,      {\Lambda} = max L({\mu}_o,{\Sigma}) / max L({\mu},{\Sigma}),    0 < {\Lambda} <= 1','Fontsize',14)
text(0.565,0.46,'{\Sigma}','Fontsize',8)
text(0.655,0.46,'{\mu},{\Sigma}','Fontsize',8)
text(0,0.2,'Since         {\Lambda} = (1+T^2/(n-1))^-^n^/^2,        rejection of H_o for large T^2 is equivalent to rejection of H_o for small {\Lambda} ','Fontsize',14)
text(0,-0.05,'and since  T^2 ~ p(n-1)/(n-p)*F_p_,_n_-_p,   the distribution of {\Lambda} is not necessary in this case to carry out the Likelihood Ratio Test ','Fontsize',14)

subplot(5,1,2:5)
Hotelling_T2 = 0:0.01:10;
n = 5;
LAMBDA_LRT = (1./(1+Hotelling_T2/(n-1))).^(n/2);
plot(Hotelling_T2,LAMBDA_LRT,'b','LineWidth',1.5),grid,xlabel('T^2','Fontsize',18),ylabel('{\Lambda}','Fontsize',24)
hold on
n = 10;
LAMBDA_LRT = (1./(1+Hotelling_T2/(n-1))).^(n/2);
plot(Hotelling_T2,LAMBDA_LRT,'c','LineWidth',1.5)
n = 30;
LAMBDA_LRT = (1./(1+Hotelling_T2/(n-1))).^(n/2);
plot(Hotelling_T2,LAMBDA_LRT,'r','LineWidth',1.5)
n = 50;
LAMBDA_LRT = (1./(1+Hotelling_T2/(n-1))).^(n/2);
plot(Hotelling_T2,LAMBDA_LRT,'g','LineWidth',1.5)
n = 100;
LAMBDA_LRT = (1./(1+Hotelling_T2/(n-1))).^(n/2);
plot(Hotelling_T2,LAMBDA_LRT,'k','LineWidth',1.5),legend('n = 5','n = 10','n = 30','n = 50','n = 100')
hold off
