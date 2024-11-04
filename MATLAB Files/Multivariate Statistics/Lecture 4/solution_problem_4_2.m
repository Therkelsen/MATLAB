clc;clear;close all
format compact
%--------------------------------------------------------------------------
% load data
%--------------------------------------------------------------------------
load("Lecture 4/dataset_problem_4_2.mat")
disp('-------------------------------------------------------------------')
disp('Bartlett test,  H0: SIGMA_female = SIGMA_male')
disp('-------------------------------------------------------------------')
%--------------------------------------------------------------------------
% descriptive statistics for Bartlett test
%--------------------------------------------------------------------------
p = 2;
n_female = 45;
n_male = 45;
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
%------------------------------------------------------------------------
% Model check for FEMALE data
%------------------------------------------------------------------------
figure
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for FEMALE observations possibly being from bivariate normal distribution (~ MVN with p = 2)','Fontsize',20)
subplot(5,2,3)
axis off
x1_aver = mean(X_female(:,1));
x2_aver = mean(X_female(:,2));
S = cov(X_female);
n = n_female;
text(0,1,'Descriptive statistics:','Fontsize',18)
text(0.5,1,['n = ' num2str(n)],'Fontsize',16)
text(0.5,0.75,['{\mu}_h_a_t = x = ' mat2str([x1_aver x2_aver],2)],'Fontsize',16)
text(0.605,0.96,'_','Fontsize',16),warning('off')
text(0.5,0.5,['{\Sigma}_h_a_t = S = ' mat2str(S,3)],'Fontsize',16)
subplot(5,2,5:2:9)
plotmatrix(X_female,'r*')
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
    d_i_sqr(i) = (X_female(i,:)-mu_hat)*inv(S)*(X_female(i,:)-mu_hat)';
end
df = 2;
z_i = chi2rnd(df,1,n);
% z2 = chi2rnd(20,1,n);
qqplot(d_i_sqr,z_i),grid,xlabel('quantiles for d_i^2','Fontsize',16),ylabel('quantiles for {\chi}_2^2 distribution','Fontsize',16),...
    title('qq-plot for d_i^2 versus {\chi}_2^2 distribution','Fontsize',16)
%------------------------------------------------------------------------
% Model check for MALE data
%------------------------------------------------------------------------
figure
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for MALE observations possibly being from bivariate normal distribution (~ MVN with p = 2)','Fontsize',20)
subplot(5,2,3)
axis off
x1_aver = mean(X_male(:,1));
x2_aver = mean(X_male(:,2));
S = cov(X_male);
n = n_male;
text(0,1,'Descriptive statistics:','Fontsize',18)
text(0.5,1,['n = ' num2str(n)],'Fontsize',16)
text(0.5,0.75,['{\mu}_h_a_t = x = ' mat2str([x1_aver x2_aver],2)],'Fontsize',16)
text(0.605,0.96,'_','Fontsize',16),warning('off')
text(0.5,0.5,['{\Sigma}_h_a_t = S = ' mat2str(S,3)],'Fontsize',16)
subplot(5,2,5:2:9)
plotmatrix(X_male,'b*')
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
    d_i_sqr(i) = (X_male(i,:)-mu_hat)*inv(S)*(X_male(i,:)-mu_hat)';
end
df = 2;
z_i = chi2rnd(df,1,n);
% z2 = chi2rnd(20,1,n);
qqplot(d_i_sqr,z_i),grid,xlabel('quantiles for d_i^2','Fontsize',16),ylabel('quantiles for {\chi}_2^2 distribution','Fontsize',16),...
    title('qq-plot for d_i^2 versus {\chi}_2^2 distribution','Fontsize',16)
%--------------------------------------------------------------------------
% further descriptive statistics for test of equal mean vectors
%--------------------------------------------------------------------------
disp('----------------------------------------------------------------------')
disp('Further descriptive statistics')
disp('----------------------------------------------------------------------')
mu_hat_female = mean(X_female)
mu_hat_male = mean(X_male)
%--------------------------------------------------------------------------
% Exact test of hypothesis   H0: mu_female = mu_male
%--------------------------------------------------------------------------
figure
subplot(4,2,8)
axis off
text(-0.2,0.8,['Exact test of hypothesis   H_0: {\mu}_f_e_m_a_l_e = {\mu}_m_a_l_e ,  {\alpha} = ' num2str(alpha)],'Fontsize',16)
T2 = (mu_hat_female-mu_hat_male)*inv(Sigma_hat_pooled*(1/n_female+1/n_male))*(mu_hat_female-mu_hat_male)'
critical_value = (p*(n_female+n_male-2)/(n_female+n_male-p-1))*finv(1-alpha,p,n_female+n_male-p-1)
text(-0.2,0.6,['T^2 =  ' num2str(T2) '  >  critical value  = (p*(n_f+n_m-2)/(n_f+n_m-p-1))*F_p_,_n_f_+_n_m_-_p_-_1({\alpha})  =  ' num2str(critical_value)],'Fontsize',14)
p_value = 1-fcdf(((n_female+n_male-p-1)/(p*(n_female+n_male-2)))*T2,p,n_female+n_male-p-1)
text(-0.2,0.45,['p-value =  ' num2str(p_value)],'Fontsize',14)
%--------------------------------------------------------------------------
% Large sample approximate test of hypothesis   H0: mu_female = mu_male
%-------------------------------------------------------------------------
text(-0.2,0.1,'Large sample approximate test of hypothesis:','Fontsize',16)
critical_value = chi2inv(1-alpha,p)
text(-0.2,-0.1,['T^2 =  ' num2str(T2) '  >  critical value  = {\chi}_2^2({\alpha})  =  ' num2str(critical_value)],'Fontsize',14)
p_value = 1-chi2cdf(T2,p)
text(-0.2,-0.25,['p-value =  ' num2str(p_value)],'Fontsize',14)
%--------------------------------------------------------------------------
% further descriptive statistics for plot of confidence region
%--------------------------------------------------------------------------
mu_diff_hat = mu_hat_female-mu_hat_male
mu_diff1_hat = mu_diff_hat(1);
mu_diff2_hat = mu_diff_hat(2);
sigma1_hat = sqrt(Sigma_hat_pooled(1,1));
sigma2_hat = sqrt(Sigma_hat_pooled(2,2));
maxvar = max(Sigma_hat_pooled(1,1),Sigma_hat_pooled(2,2));
%--------------------------------------------------------------------------
% Confidence region for mu_female - mu_male (not for observations !)
%--------------------------------------------------------------------------
subplot(4,2,2:2:6)
y1 = mu_diff1_hat-sqrt(Sigma_hat_pooled(1,1)):sqrt(Sigma_hat_pooled(1,1))/50:mu_diff1_hat+sqrt(Sigma_hat_pooled(1,1)); 
y2 = mu_diff2_hat-sqrt(Sigma_hat_pooled(2,2)):sqrt(Sigma_hat_pooled(2,2))/50:mu_diff2_hat+sqrt(Sigma_hat_pooled(2,2));
[Y1,Y2] = meshgrid(y1,y2);
F = mvnpdf([Y1(:) Y2(:)],mu_diff_hat,Sigma_hat_pooled);
F = reshape(F,length(y2),length(y1));
alpha = 0.05;
alpha_contour = [0.05 0.05];
%c = (p*(n-1)/(n*(n-p)))*finv(1-alpha,p,n-p);
c = (p*(n_female+n_male-2)/(n_female+n_male-p-1))*(1/n_female+1/n_male)*finv(1-alpha_contour,p,n_female+n_male-p-1);
C = (1/(2*pi*sqrt(det(Sigma_hat_pooled))))*exp(-c/2);
%contour(y1,y2,F,C,'Color','k','LineWidth',3),grid,axis equal,xlabel('{\mu}_f_e_m_a_l_e_,_1 - {\mu}_m_a_l_e_,_1','Fontsize',18),ylabel('{\mu}_f_e_m_a_l_e_,_2 - {\mu}_m_a_l_e_,_2','Fontsize',18),...
        title({'95% Confidence region and intervals for ({\mu}_f_e_m_a_l_e_,_1-{\mu}_m_a_l_e_,_1, {\mu}_f_e_m_a_l_e_,_2-{\mu}_m_a_l_e_,_2 based on observations:';'CR (black), simult.CI (red), marg.CI (blue), Bonf.CI (green)'},'Fontsize',14)
plot2d_CR_for_difference_in_mu_ellipsis (mu_diff_hat', Sigma_hat_pooled, alpha, n_female, n_male)
hold on
%--------------------------------------------------------------------------
% Simultaneous confidence intervals for mu_female - mu_male (not for observations !)
%--------------------------------------------------------------------------
disp('----------------------------------------------------------------------')
disp('Confidence intervals for mu_female - mu_male')
disp('----------------------------------------------------------------------')
mu_diff1_sim_CI = mu_diff1_hat + [-1 1]*sqrt((p*(n_female+n_male-2)/(n_female+n_male-p-1))*finv(1-alpha,p,n_female+n_male-p-1))*sqrt(Sigma_hat_pooled(1,1)*(1/n_female+1/n_male))
line([mu_diff1_sim_CI(1) mu_diff1_sim_CI(1)],[-12 10],'LineWidth',2,'Color','r')
line([mu_diff1_sim_CI(2) mu_diff1_sim_CI(2)],[-12 10],'LineWidth',2,'Color','r')
mu_diff2_sim_CI = mu_diff2_hat + [-1 1]*sqrt((p*(n_female+n_male-2)/(n_female+n_male-p-1))*finv(1-alpha,p,n_female+n_male-p-1))*sqrt(Sigma_hat_pooled(2,2)*(1/n_female+1/n_male))
line([-1 15],[mu_diff2_sim_CI(1) mu_diff2_sim_CI(1)],'LineWidth',2,'Color','r')
line([-1 15],[mu_diff2_sim_CI(2) mu_diff2_sim_CI(2)],'LineWidth',2,'Color','r')
%--------------------------------------------------------------------------
% Marginal confidence intervals for mu_female - mu_male (not for observations !)
%--------------------------------------------------------------------------
mu_diff1_marg_CI = mu_diff1_hat + [-1 1]*tinv(1-alpha/2,n_female+n_male-2)*sqrt(Sigma_hat_pooled(1,1)*(1/n_female+1/n_male))
line([mu_diff1_marg_CI(1) mu_diff1_marg_CI(1)],[-12 10],'LineWidth',2,'Color','b')
line([mu_diff1_marg_CI(2) mu_diff1_marg_CI(2)],[-12 10],'LineWidth',2,'Color','b')
mu_diff2_marg_CI = mu_diff2_hat + [-1 1]*tinv(1-alpha/2,n_female+n_male-2)*sqrt(Sigma_hat_pooled(2,2)*(1/n_female+1/n_male))
line([-1 15],[mu_diff2_marg_CI(1) mu_diff2_marg_CI(1)],'LineWidth',2,'Color','b')
line([-1 15],[mu_diff2_marg_CI(2) mu_diff2_marg_CI(2)],'LineWidth',2,'Color','b')
%--------------------------------------------------------------------------
% Bonferroni simultaneous confidence intervals for mu_female - mu_male (not for observations !)
%--------------------------------------------------------------------------
mu_diff1_bonf_CI = mu_diff1_hat + [-1 1]*tinv(1-alpha/(2*p),n_female+n_male-2)*sqrt(Sigma_hat_pooled(1,1)*(1/n_female+1/n_male))
line([mu_diff1_bonf_CI(1) mu_diff1_bonf_CI(1)],[-12 10],'LineWidth',2,'Color','g')
line([mu_diff1_bonf_CI(2) mu_diff1_bonf_CI(2)],[-12 10],'LineWidth',2,'Color','g')
mu_diff2_bonf_CI = mu_diff2_hat + [-1 1]*tinv(1-alpha/(2*p),n_female+n_male-2)*sqrt(Sigma_hat_pooled(2,2)*(1/n_female+1/n_male))
line([-1 15],[mu_diff2_bonf_CI(1) mu_diff2_bonf_CI(1)],'LineWidth',2,'Color','g')
line([-1 15],[mu_diff2_bonf_CI(2) mu_diff2_bonf_CI(2)],'LineWidth',2,'Color','g')
hold off
disp('----------------------------------------------------------------------')
%--------------------------------------------------------------------------
% Text for figure
%--------------------------------------------------------------------------
subplot(4,2,1)
axis off
text(-0.3,0.4,'Confidence region (CR) and confidence intervals (CI) for','Fontsize',20)
text(-0.3,0.1,'{\mu}_f_e_m_a_l_e -  {\mu}_m_a_l_e  for two assumed bivariate normal distributions','Fontsize',20)
text(-0.3,-0.2,'X_f_e_m_a_l_e~ N_2({\mu}_f_e_m_a_l_e,{\Sigma}_f_e_m_a_l_e)  and  X_m_a_l_e~ N_2({\mu}_m_a_l_e,{\Sigma}_m_a_l_e)','Fontsize',20)
text(-0.3,-0.95,'See command window output for','Fontsize',20)
text(-0.2,-1.3,'- Bartlett test of  {\Sigma}_f_e_m_a_l_e= {\Sigma}_m_a_l_e','Fontsize',18)
text(-0.2,-1.6,'- Numerical descriptive statistics','Fontsize',18)
text(-0.2,-1.9,'- Numerical simultaneous confidence intervals for {\mu}_f_e_m_a_l_e -  {\mu}_m_a_l_e ','Fontsize',18)
text(-0.2,-2.2,'- Numerical marginal confidence intervals for {\mu}_f_e_m_a_l_e -  {\mu}_m_a_l_e ','Fontsize',18)
text(-0.2,-2.5,'- Numerical Bonferroni confidence intervals for {\mu}_f_e_m_a_l_e -  {\mu}_m_a_l_e ','Fontsize',18)
%--------------------------------------------------------