clc;clear;close all
format compact
utils = Utils;
%--------------------------------------------------------------------------
% load data
%--------------------------------------------------------------------------
load("Lecture 3/dataset_problem_3_1.mat");
%--------------------------------------------------------------------------
% descriptive statistics
%--------------------------------------------------------------------------
p = 2;
n = 30;
mu_hat = mean(obs);
mu1_hat = mu_hat(1);
mu2_hat = mu_hat(2);
Sigma_hat = cov(obs);
sigma1_hat = sqrt(Sigma_hat(1,1));
sigma2_hat = sqrt(Sigma_hat(2,2));
maxvar = max(Sigma_hat(1,1),Sigma_hat(2,2));
subplot(4,2,3)
axis off
%--------------------------------------------------------------------------
%  plot observations and prediction ellipses for observations
%--------------------------------------------------------------------------
x1 = mu1_hat-3*sigma1_hat:6*sigma1_hat/50:mu1_hat+3*sigma1_hat; 
x2 = mu2_hat-3*sigma2_hat:6*sigma2_hat/50:mu2_hat+3*sigma2_hat;
[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],mu_hat,Sigma_hat);
F = reshape(F,length(x2),length(x1));
subplot(4,2,[5 7])
alpha = [0.5 0.25 0.05];
plot(obs(:,1),obs(:,2),'+'),grid,xlabel('x_1','Fontsize',16); ylabel('x_2','Fontsize',16),...
    title({['n = ' num2str(n) ' observations scatter plot for (X_1, X_2) and'];['(1-{\alpha})*100% prediction ellipses for  {\alpha}  =  [' num2str(alpha(3)) '  ' num2str(alpha(2)) '  ' num2str(alpha(1)) ']']},'Fontsize',14)
hold on
c = chi2inv(1-alpha,p);
C = (1/(2*pi*sqrt(det(Sigma_hat)))).*exp(-c/2);
contour(x1,x2,F,C,'Color','k','LineWidth',2),axis square
text(-2000,1000,'Descriptive statistics:','Fontsize',14)
text(-2000,900,['{\mu}_h_a_t = x_m_e_a_n = [' num2str(mu1_hat,3) '  ' num2str(mu2_hat,3) ']'],'Fontsize',12)
text(-2000,800,['{\Sigma}_h_a_t = S = [' num2str(Sigma_hat(1,1),3) '  ' num2str(Sigma_hat(1,2),3) ';'],'Fontsize',12)
text(-2000,750,['                  ' num2str(Sigma_hat(2,1),3) '  ' num2str(Sigma_hat(2,2),3) ']'],'Fontsize',12)
hold off
%--------------------------------------------------------------------------
% Confidence region for mu (not for observations !)
%--------------------------------------------------------------------------
subplot(4,2,2:2:6)
y1 = mu1_hat-sqrt(Sigma_hat(1,1)):sqrt(Sigma_hat(1,1))/50:mu1_hat+sqrt(Sigma_hat(1,1)); 
y2 = mu2_hat-sqrt(Sigma_hat(2,2)):sqrt(Sigma_hat(2,2))/50:mu2_hat+sqrt(Sigma_hat(2,2));
[Y1,Y2] = meshgrid(y1,y2);
F = mvnpdf([Y1(:) Y2(:)],mu_hat,Sigma_hat);
F = reshape(F,length(y2),length(y1));
alpha = 0.05;
alpha_contour = [0.05 0.05];
% c = (p*(n-1)/(n*(n-p)))*finv(1-alpha_contour,p,n-p);
% C = (1/(2*pi*sqrt(det(Sigma_hat))))*exp(-c/2);
% contour(y1,y2,F,C,'Color','k','LineWidth',3),grid,xlabel('{\mu}_1','Fontsize',18),ylabel('{\mu}_2','Fontsize',18),xlim([1650 2050]),ylim([650 1000]),...
        title({'95% Confidence region and intervals for ({\mu}_1, {\mu}_2) based on observations:';'CR (black), simult.CI (red), marg.CI (blue), Bonf.CI (green)'},'Fontsize',14)
utils.plot2d_CR_for_mu_ellipsis (mu_hat',Sigma_hat,alpha,n)
title({'95% Confidence region and intervals for ({\mu}_1, {\mu}_2) based on observations:';'CR (black), simult.CI (red), marg.CI (blue), Bonf.CI (green)'},'Fontsize',14)
hold on
[V, lambda] = eig(Sigma_hat);
ellipse_axes_half_lengths_1 = sqrt(lambda(1,1))*sqrt((p*(n-1)/(n*(n-p)))*finv(1-alpha,p,n-p));
ellipse_axes_half_lengths_2 = sqrt(lambda(2,2))*sqrt((p*(n-1)/(n*(n-p)))*finv(1-alpha,p,n-p));
line([mu1_hat  mu1_hat-ellipse_axes_half_lengths_1*V(1,1)],[mu2_hat mu2_hat-ellipse_axes_half_lengths_1*V(2,1)],'LineWidth',2,'Color','k')
line([mu1_hat  mu1_hat+ellipse_axes_half_lengths_2*V(1,2)],[mu2_hat mu2_hat+ellipse_axes_half_lengths_2*V(2,2)],'LineWidth',2,'Color','k')
%--------------------------------------------------------------------------
% Simultaneous confidence intervals for mu (not for observations !)
%--------------------------------------------------------------------------
mu1_sim_CI = [mu1_hat-sqrt((p*(n-1)/(n-p))*finv(1-alpha,p,n-p))*sqrt(Sigma_hat(1,1)/n) mu1_hat+sqrt((p*(n-1)/(n-p))*finv(1-alpha,p,n-p))*sqrt(Sigma_hat(1,1)/n)]
line([mu1_sim_CI(1) mu1_sim_CI(1)],[650 1000],'LineWidth',2,'Color','r')
line([mu1_sim_CI(2) mu1_sim_CI(2)],[650 1000],'LineWidth',2,'Color','r')
mu2_sim_CI = [mu2_hat-sqrt((p*(n-1)/(n-p))*finv(1-alpha,p,n-p))*sqrt(Sigma_hat(2,2)/n) mu2_hat+sqrt((p*(n-1)/(n-p))*finv(1-alpha,p,n-p))*sqrt(Sigma_hat(2,2)/n)]
line([1650 2050],[mu2_sim_CI(1) mu2_sim_CI(1)],'LineWidth',2,'Color','r')
line([1650 2050],[mu2_sim_CI(2) mu2_sim_CI(2)],'LineWidth',2,'Color','r')
%--------------------------------------------------------------------------
% Marginal confidence intervals for mu (not for observations !)
%--------------------------------------------------------------------------
mu1_marg_CI = [mu1_hat-tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(1,1)/n) mu1_hat+tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(1,1)/n)]
line([mu1_marg_CI(1) mu1_marg_CI(1)],[650 1000],'LineWidth',2,'Color','b')
line([mu1_marg_CI(2) mu1_marg_CI(2)],[650 1000],'LineWidth',2,'Color','b')
mu2_marg_CI = [mu2_hat-tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(2,2)/n) mu2_hat+tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(2,2)/n)]
line([1650 2050],[mu2_marg_CI(1) mu2_marg_CI(1)],'LineWidth',2,'Color','b')
line([1650 2050],[mu2_marg_CI(2) mu2_marg_CI(2)],'LineWidth',2,'Color','b')
%--------------------------------------------------------------------------
% Bonferroni simultaneous confidence intervals for mu (not for observations !)
%--------------------------------------------------------------------------
mu1_bonf_CI = [mu1_hat-tinv(1-alpha/(2*p),n-1)*sqrt(Sigma_hat(1,1)/n) mu1_hat+tinv(1-alpha/(2*p),n-1)*sqrt(Sigma_hat(1,1)/n)]
line([mu1_bonf_CI(1) mu1_bonf_CI(1)],[650 1000],'LineWidth',2,'Color','g')
line([mu1_bonf_CI(2) mu1_bonf_CI(2)],[650 1000],'LineWidth',2,'Color','g')
mu2_bonf_CI = [mu2_hat-tinv(1-alpha/(2*p),n-1)*sqrt(Sigma_hat(2,2)/n) mu2_hat+tinv(1-alpha/(2*p),n-1)*sqrt(Sigma_hat(2,2)/n)]
line([1650 2050],[mu2_bonf_CI(1) mu2_bonf_CI(1)],'LineWidth',2,'Color','g')
line([1650 2050],[mu2_bonf_CI(2) mu2_bonf_CI(2)],'LineWidth',2,'Color','g')
text(0.999*mu1_hat,mu2_hat,'o','Fontsize',14)
text(mu1_hat-0.075,mu2_hat+0.03,'({\mu}_1_,_h_a_t, {\mu}_2_,_h_a_t)','Fontsize',14)
hold off
%--------------------------------------------------------------------------
% Exact test of hypothesis   H0: mu = mu0
%--------------------------------------------------------------------------
subplot(4,2,8)
axis off
mu0 = [1750 950];
text(0,0.8,['Exact test of hypothesis   H_0: {\mu} = {\mu}_0 = ' mat2str(mu0) ',  {\alpha} = ' num2str(alpha)],'Fontsize',16)
T2 = (mu_hat-mu0)*inv(Sigma_hat/n)*(mu_hat-mu0)'
critical_value = (p*(n-1)/(n-p))*finv(1-alpha,p,n-p)
text(0,0.6,['T^2 =  ' num2str(T2) '  >  critical value  = p(n-1)/(n-p)*F_p_,_n_-_p({\alpha})  =  ' num2str(critical_value)],'Fontsize',14)
p_value = 1-fcdf(((n-p)/(p*(n-1)))*T2,p,n-p);
text(0,0.45,['p-value =  ' num2str(p_value)],'Fontsize',14)
%--------------------------------------------------------------------------
% Large sample approximate test of hypothesis   H0: mu = mu0
%-------------------------------------------------------------------------
text(0,0.1,'Large sample approximate test of hypothesis:','Fontsize',16)
critical_value = chi2inv(1-alpha,p)
text(0,-0.1,['T^2 =  ' num2str(T2) '  >  critical value  = {\chi}_2^2({\alpha})  =  ' num2str(critical_value)],'Fontsize',14)
p_value = 1-chi2cdf(T2,p);
text(0,-0.25,['p-value =  ' num2str(p_value)],'Fontsize',14)
%--------------------------------------------------------------------------
% LRT approximate test of hypothesis   H0: mu = mu0
%-------------------------------------------------------------------------
Sigma_hat_0 = (obs-mu0)'*(obs-mu0)/(n-1)
LAMBDA_LR = ((det(Sigma_hat)/det(Sigma_hat_0)))^(n/2)
LAMBDA_LR_from_T2 = (1 + T2/(n-1))^(-n/2)
LRT_TS = -2*log(LAMBDA_LR)
%--------------------------------------------------------------------------
% Text for figure
%--------------------------------------------------------------------------
subplot(4,2,1)
axis off
text(-0.3,1,'Confidence region (CR) and confidence intervals (CI) for {\mu} = ({\mu}_1,{\mu}_2)^T','Fontsize',20)
text(-0.3,0.75,'of an assumed bivariate normal, (X_1,X_2) ~ N_2({\mu},{\Sigma})   ({\mu},{\Sigma} unknown)','Fontsize',20)
text(-0.3,0.5,'Level of confidence:  100(1-{\alpha}) %','Fontsize',16)
text(-0.3,0.2,'CR (two-dimens.):     (x-{\mu})^T*S^-^1*(x-{\mu})  {\leq}  p(n-1)/n(n-p) * F_p_,_n_-_p({\alpha})','Fontsize',16)
text(-0.095,0.37,'                 _              _','Fontsize',16)
text(-0.3,-0.05,'Simultaneous CI''s:    ( x_i   {\pm}  sqrt [ p(n-1)/(n-p) * F_p_,_n_-_p({\alpha}) ] * sqrt [ s_i_i/n ] ),   i=1,2','Fontsize',16)
text(-0.093,0.13,'                 _','Fontsize',16)
text(-0.3,-0.3,'Marginal CI''s:            ( x_i   {\pm}  t_n_-_1({\alpha}/2) * sqrt [ s_i_i/n ] ),   i=1,2','Fontsize',16)
text(-0.093,-0.12,'                 _','Fontsize',16)
text(-0.3,-0.55,'Bonferroni CI''s:         ( x_i   {\pm}  t_n_-_1({\alpha}/2p) * sqrt [ s_i_i/n ] ),   i=1,2','Fontsize',16)
text(-0.093,-0.37,'                 _','Fontsize',16)
%------------------------------------------------------------------------------------------
% Model check
%------------------------------------------------------------------------
figure
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for observations possibly being from bivariate normal distribution (~ MVN with p = 2)','Fontsize',20)
subplot(5,2,3)
axis off
text(0,1.2,'Descriptive statistics:','Fontsize',18)
text(0.5,1.2,['n = ' num2str(n)],'Fontsize',16)
text(0.5,0.85,['{\mu}_h_a_t = x = ' mat2str([mu1_hat mu2_hat],3)],'Fontsize',16)
text(0.605,1.07,'_','Fontsize',16),warning('off')
text(0.5,0.5,['{\Sigma}_h_a_t = S = ' mat2str(Sigma_hat,3)],'Fontsize',16)
subplot(5,2,5:2:9)
plotmatrix(obs,'r*')
xlabel('x_1                                          x_2','Fontsize',16)
ylabel('x_2                           x_1','Fontsize',16)
title(['Scattermatrix of observations'],'Fontsize',16)
subplot(5,2,4)
axis off
text(0,1.2,'Sample Mahalanobis distances, i=1..n :','Fontsize',18)
text(0,0.85,'d_i^2 = (x_i - x)^T S^-^1 (x_i - x) ~ {\chi}_2^2   (approximately)','Fontsize',16)
text(0.158,1.05,'_','Fontsize',16),warning('off')
text(0.358,1.05,'_','Fontsize',16),warning('off')
subplot(5,2,6:2:10)
mu_hat = [mu1_hat mu2_hat];
d_i_sqr = zeros(1,n);
for i = 1:n
    d_i_sqr(i) = (obs(i,:)-mu_hat)*inv(Sigma_hat)*(obs(i,:)-mu_hat)';
end
df = 2;
z_i = chi2rnd(df,1,n);
qqplot(d_i_sqr,z_i),grid,xlabel('quantiles for d_i^2','Fontsize',16),ylabel('quantiles for {\chi}_2^2 distribution','Fontsize',16),...
    title('qq-plot for d_i^2 versus {\chi}_2^2 distribution','Fontsize',16)
%------------------------------------------------------------------------