clc;clear;close all
format compact 

subplot(5,2,1:2)
axis off
text(0,1.35,'LARGE SAMPLE approximate hypothesis testing for multivariate (not necessarily normal) population mean vector,  H_o: {\mu} = {\mu}_o','Fontsize',16)
text(0,1,'For n >> p, the Hotelling test-statistic is due to CLT APPROXIMATELY Chi2-distributed,                  T^2 = (x - {\mu}_o)^T (S/n)^-^1 (x - {\mu}_o) ~ {\chi}^2_p','Fontsize',14)
text(0.82,1.19,'_','Fontsize',16),warning('off')
text(0.71,1.19,'_','Fontsize',16),warning('off')
text(0,0.82,'even without assuming a multivariate normal distribution for X','Fontsize',14) 
text(0,0.5,'If, in fact, X is multivariate NORMAL distributed, the EXACT distribution of the test statistic is           T^2 = (x - {\mu}_o)^T (S/n)^-^1 (x - {\mu}_o) ~ p(n-1)/(n-p)*F_p_,_n_-_p ','Fontsize',14)
text(0.82,0.69,'_','Fontsize',16),warning('off')
text(0.71,0.69,'_','Fontsize',16),warning('off')
text(0,0.2,'Here the approximate large-sample equivalence of these two test-distributions is assessed by comparing quantiles (critical values), {\chi}^2_p({\alpha}) and (p(n-1)/(n-p))F_p_,_n_-_p({\alpha})','Fontsize',14)

%--------------------------------------------------------------------------
% distributions
%--------------------------------------------------------------------------
p = 5;
df_chi2 = p;
n = 100;
df_F_1 = p;
df_F_2 = n-p;

x = 0:0.01:16;
f_chi2 = chi2pdf(x,df_chi2);
f_F = fpdf(x,df_F_1,df_F_2);

alfa = 0.1;
quantile_chi2 = chi2inv(1-alfa,df_chi2)
quantile_F = finv(1-alfa,df_F_1,df_F_2)
scaled_quantile_F = (p*(n-1)/(n-p))*quantile_F

subplot(5,2,3:2:9)
plot(x,f_chi2,'b','LineWidth',2),grid,xlim([0 16]),xlabel('x','Fontsize',16),ylabel('pdf','Fontsize',16)
hold on
line([quantile_chi2 quantile_chi2],[0 0.2],'Color','b','LineWidth',3)
plot(x,f_F,'r','LineWidth',2)
line([quantile_F quantile_F],[0 0.4],'Color','r','LineWidth',3)
line([scaled_quantile_F scaled_quantile_F],[0 0.3],'Color','r','LineWidth',3)
hold off
text(4.2,0.83,'distributions and quantiles','Fontsize',16)
text(1.1,0.72,'pdf F_p_,_n_-_p(x)','Fontsize',14,'color','r')
text(1.5,0.47,'quantile','Fontsize',14,'color','r')
text(1.4,0.43,'F_p_,_n_-_p({\alpha})=1.91','Fontsize',14,'color','r')
text(8.2,0.36,'scaled quantile','Fontsize',14,'color','r')
text(7,0.32,'(p(n-1)/(n-p))F_p_,_n_-_p({\alpha})=9.95','Fontsize',14,'color','r')
text(3.1,0.18,'pdf {\chi}^2_p(x)','Fontsize',14,'color','b')
text(7.5,0.26,'quantile','Fontsize',14,'color','b')
text(7,0.22,'{\chi}^2_p({\alpha})=9.24','Fontsize',14,'color','b')
text(8.1,0.72,'Example parameters:','Fontsize',14,'color','k')
text(8.1,0.68,'p = 5','Fontsize',14,'color','k')
text(8.1,0.65,'n = 100 >> p','Fontsize',14,'color','k')
text(8.1,0.62,'{\alpha} = 0.1','Fontsize',14,'color','k')

%--------------------------------------------------------------------------
% quantile ratios
%--------------------------------------------------------------------------
subplot(5,2,4:2:10)

alfa1 = 0.1;
p = 5;
df_chi2 = p;
quantile_chi2 = chi2inv(1-alfa1,df_chi2);

n = 30:1000;
df_F_1 = p;
df_F_2 = n-p;
scaled_quantile_F = (p*(n-1)./(n-p)).*finv(1-alfa1,df_F_1,df_F_2);

quantiles_ratio_alfa1 = quantile_chi2./scaled_quantile_F;
plot(n,quantiles_ratio_alfa1,'b','LineWidth',2),grid,xlabel('sample size, n','Fontsize',16),ylabel('quantile ratio','Fontsize',16)
text(110,1.015,'quantile ratio,  {\chi}^2_p({\alpha}) / [(p(n-1)/(n-p))F_p_,_n_-_p({\alpha})]','Fontsize',16)
text(610,0.88,'Example parameter:','Fontsize',14,'color','k')
text(610,0.865,'p = 5','Fontsize',14,'color','k')
hold on

alfa2 = 0.05;
quantile_chi2 = chi2inv(1-alfa2,df_chi2);
scaled_quantile_F = (p*(n-1)./(n-p)).*finv(1-alfa2,df_F_1,df_F_2);
quantiles_ratio_alfa2 = quantile_chi2./scaled_quantile_F;
plot(n,quantiles_ratio_alfa2,'r','LineWidth',2)

alfa3 = 0.01;
quantile_chi2 = chi2inv(1-alfa3,df_chi2);
scaled_quantile_F = (p*(n-1)./(n-p)).*finv(1-alfa3,df_F_1,df_F_2);
quantiles_ratio_alfa3 = quantile_chi2./scaled_quantile_F;
plot(n,quantiles_ratio_alfa3,'g','LineWidth',2)

alfa4 = 0.001;
quantile_chi2 = chi2inv(1-alfa4,df_chi2);
scaled_quantile_F = (p*(n-1)./(n-p)).*finv(1-alfa4,df_F_1,df_F_2);
quantiles_ratio_alfa4 = quantile_chi2./scaled_quantile_F;
plot(n,quantiles_ratio_alfa4,'k','LineWidth',2),legend(['{\alpha} = ' num2str(alfa1)],['{\alpha} = ' num2str(alfa2)],['{\alpha} = ' num2str(alfa3)],['{\alpha} = ' num2str(alfa3)],'Location','East')
hold off

