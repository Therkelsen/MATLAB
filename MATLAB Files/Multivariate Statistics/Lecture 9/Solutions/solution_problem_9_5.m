clc,clear,close all
format compact
%----------------------------------------
% data 
%----------------------------------------
load('C:\Users\cv\OneDrive - Syddansk Universitet\M-drive moved by SDU-IT\Statistik\MULTISTAT\Lektion 09  Non-parametric statistics\Course material\dataset_problem_9_5.mat')
disp('-----------------------------------')
disp('Goodness Of Fit test')
disp('-----------------------------------')
n = length(X)
%----------------------------------------
% hypothesized distribution 
%----------------------------------------
disp('---------------------------------------------')
disp('hypothesized distribution')
disp('---------------------------------------------')
disp('Exponential distribution')
lambda_hat = 1/mean(X)
%----------------------------------------
% visualization of data vs. model
%----------------------------------------
figure(1)
H0_data = exprnd(1/lambda_hat,1,1e6);
qqplot(X,H0_data),title('QQplot data vs. model'),xlabel('quantiles for estimated model'),ylabel('data quantiles')
%------------------------------------
figure(2)
hist_bins = 10;
histogram(X,hist_bins,'Normalization','pdf');
hold on
x = 0:0.01:5/lambda_hat;
X_0_pdf = exppdf(x,1/lambda_hat);       
plot(x,X_0_pdf,'r','Linewidth',2),title('data vs. model pdf')   
legend(' data histogram',' pdf for estimated model','Location','NorthEast')
hold off
%--------------------------------------
figure(3)
H = cdfplot(X);
H.LineWidth = 1.5;
hold on
X_0_cdf = expcdf(x,1/lambda_hat);       
plot(x,X_0_cdf,'r','Linewidth',1.5),title('data vs. model cdf')  
legend(' data cdf',' cdf for estimated model','Location','East')
hold off
%----------------------------------------
% test statistic T (Pearson's GOF X^2)
%----------------------------------------
disp('---------------------------------------------')
disp('test statistic T (Pearson''s GOF X^2)')
disp('---------------------------------------------')
Nbins = 10
bin_limits = expinv(0:1/Nbins:1,1/lambda_hat)
Ei = n/Nbins
Oi = zeros(1,Nbins);
for i = 1:Nbins
    Oi(i) = length(find(X > bin_limits(i) & X < bin_limits(i+1)));
end
Oi
T = sum((Oi-Ei).^2)/Ei
%-------------------------------------
figure(4)
for i = 1:Nbins-1
    line([bin_limits(i) bin_limits(i)],[0 Oi(i)],'Linewidth',3,'Color','b')
    line([bin_limits(i+1) bin_limits(i+1)],[0 Oi(i)],'Linewidth',3,'Color','b')
    line([bin_limits(i) bin_limits(i+1)],[Oi(i) Oi(i)],'Linewidth',3,'Color','b')
end
line([bin_limits(Nbins) 1.25*bin_limits(Nbins)],[Oi(Nbins) Oi(Nbins)],'Linewidth',3,'Color','b')
line([0 1.25*bin_limits(Nbins)],[Ei Ei],'LineStyle','--','Linewidth',4,'Color','r')    
title('Binning for Pearson''s {\chi}^2 test:   n = 100,   Nbins = 10,   Ei (RED),   Oi(BLUE)','FontSize',16)
xlim([0 1.25*bin_limits(Nbins)])
%---------------------------------------------
% approximate test with T (Pearson''s GOF X^2)
%---------------------------------------------
disp('-----------------------------------------------------')
disp('approximate test with T (Pearson''s GOF X^2)')
disp('-----------------------------------------------------')
alpha = 0.05
Npar = 1
df = Nbins-Npar-1
critical_value = chi2inv(1-alpha,df)
reject_H0 = T > critical_value
p_value = 1-chi2cdf(T,df)
%----------------------------------------
% approximate test with built-in command
%----------------------------------------
disp('---------------------------------------------')
disp('approximate test with built-in command (Pearson''s GOF X^2)')
disp('---------------------------------------------')
prob_dist_0 = fitdist(X','Exponential')
[reject_H0,p_value,stats] = chi2gof(X,'CDF',prob_dist_0,'Edges',bin_limits)
disp('---------------------------------------------')
%----------------------------------------
