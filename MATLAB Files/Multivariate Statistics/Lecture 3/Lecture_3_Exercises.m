clc; clear; close all;
format compact

utils = Utils;

disp("Problem 3.1")

%--------------------------------------------------------------------------
% load data
%--------------------------------------------------------------------------
data = load("Lecture 3/dataset_problem_3_1.mat").obs
%--------------------------------------------------------------------------
% descriptive statistics
%--------------------------------------------------------------------------
% % Calculate descriptive statustics (Mean vector, Covariance matrix, range vector)
[mu_hat, Sigma_hat, range] = utils.calculate_descriptive_statistics(data, true, false);
% % Column 1 and 2 of the mean vector
mu1_hat = mu_hat(1); mu2_hat = mu_hat(2);
% % Standard deviations of diagonals of covariance matrix (std = sqrt(cov))
sigma1_hat = sqrt(Sigma_hat(1,1)); sigma2_hat = sqrt(Sigma_hat(2,2));
% % The biggest covariance element of the diagonals
maxvar = max(Sigma_hat(1,1),Sigma_hat(2,2));

% % Plot the observations using a scatter plot
% % n rows p columns, each p is a variable
[n, p] = size(data)

%--------------------------------------------------------------------------
%  plot observations and !!PREDICTION!! ellipses for observations
%--------------------------------------------------------------------------
% % These are vectors representing a range of values around the means of
% % the first and second variables, respectively. They cover the range from
% % 3 standard deviations below the mean to 3 standard deviations above the mean.
% 
% % For example,
% % x1 = mu1_hat-3*sigma1_hat:6*sigma1_hat/50:mu1_hat+3*sigma1_hat
% % defines a range of 50 evenly spaced points between mu1_hat - 3 * sigma1_hat
% % and mu1_hat + 3 * sigma1_hat, where 6*sigma1_hat/50 controls the step size.
% % Similarly, x2 is for the second variable.
x1 = mu1_hat-3*sigma1_hat : 6*sigma1_hat/50 : mu1_hat+3*sigma1_hat; 
x2 = mu2_hat-3*sigma2_hat : 6*sigma2_hat/50 : mu2_hat+3*sigma2_hat;

[X1,X2] = meshgrid(x1,x2);
% % MVN PDF is a multivariate normal probability density function
% % I.e. multivariate gaussian distribution
% % This is the joint distribution between the two variables.
F = mvnpdf([X1(:) X2(:)],mu_hat,Sigma_hat);
F = reshape(F,length(x2),length(x1));
figure(4)

% % Confidence levels for ellipses, i.e. 1 minus alpha.
% % Example, alpha = 0.05 -> confidence = 1 - 0.05 = 95%.
alpha = [0.5 0.25 0.05];
% % Scatter plot
plot(data(:,1),data(:,2),'+'),grid,xlabel('x_1','Fontsize',16); ylabel('x_2','Fontsize',16),...
    title({['n = ' num2str(n) ' observations scatter plot for (X_1, X_2) and'];['(1-{\alpha})*100% prediction ellipses for  {\alpha}  =  [' num2str(alpha(3)) '  ' num2str(alpha(2)) '  ' num2str(alpha(1)) ']']},'Fontsize',14)
hold on

% % calculates the inverse of the chi-squared cumulative distribution
% % function (CDF) with p degrees of freedom, which is used to determine
% % the confidence region boundaries for the multivariate normal distribution.
% % For each value in alpha, this calculates the corresponding chi-squared value,
% % c, which will be used to draw the prediction ellipses.
c = chi2inv(1-alpha,p);

% % This calculates a constant scaling factor C for the contour lines of the ellipses.
% % The formula is based on the probability density function of the multivariate normal
% % distribution and the chi-squared critical values c.
C = (1/(2*pi*sqrt(det(Sigma_hat)))).*exp(-c/2);

% % This draws the prediction ellipses on the plot. The ellipses represent
% % the 1−α confidence regions for the multivariate normal distribution.
% % x1 and x2 define the grid points, F contains the probability density values
% % and C defines the contour levels.
contour(x1,x2,F,C,'Color','k','LineWidth',2),axis square
text(-2000,1000,'Descriptive statistics:','Fontsize',14)
text(-2000,900,['{\mu}_h_a_t = x_m_e_a_n = [' num2str(mu1_hat,3) '  ' num2str(mu2_hat,3) ']'],'Fontsize',12)
text(-2000,800,['{\Sigma}_h_a_t = S = [' num2str(Sigma_hat(1,1),3) '  ' num2str(Sigma_hat(1,2),3) ';'],'Fontsize',12)
text(-2000,750,['                  ' num2str(Sigma_hat(2,1),3) '  ' num2str(Sigma_hat(2,2),3) ']'],'Fontsize',12)
hold off

%--------------------------------------------------------------------------
% Confidence region for mu (not for observations !)
%--------------------------------------------------------------------------
% % This is a vector representing a range of values around the mean of the
% % first variable (X1). It covers the range from one standard deviation below
% % the mean (𝜇1) to one standard deviation above the mean.
% % The step size is controlled by sqrt(Sigma_hat(1,1))/50,
% % creating 50 evenly spaced points. Similarly for Y2, just using X2 and 𝜇2.
y1 = mu1_hat-sqrt(Sigma_hat(1,1)) : sqrt(Sigma_hat(1,1))/50 : mu1_hat+sqrt(Sigma_hat(1,1)); 
y2 = mu2_hat-sqrt(Sigma_hat(2,2)) : sqrt(Sigma_hat(2,2))/50 : mu2_hat+sqrt(Sigma_hat(2,2));

% % MVN PDF is a multivariate normal probability density function
% % I.e. multivariate gaussian distribution
% % This is the joint distribution between the two variables.
[Y1,Y2] = meshgrid(y1,y2);
F = mvnpdf([Y1(:) Y2(:)],mu_hat,Sigma_hat);
F = reshape(F,length(y2),length(y1));

% % Confidence levels for ellipses, i.e. 1 minus alpha.
% % Example, alpha = 0.05 -> confidence = 1 - 0.05 = 95%.
alpha = 0.05;
alpha_contour = [0.05 0.05];

% % Draw the ellipses using Claus' handy little tool
figure(5)
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
% This is wide / pessimistic
%--------------------------------------------------------------------------
% % These are the confidence intervals for EACH mu!
% % So, upper and lower bound, for each mu!
% % Same thing in the two next sections.

% % There will be uncertainties in the confidence intervals, this is
% % because compounding 0.95 results in a worse result, for example 0.95^5 ≈0.77.
mu1_sim_CI = [mu1_hat-sqrt((p*(n-1)/(n-p))*finv(1-alpha,p,n-p))*sqrt(Sigma_hat(1,1)/n) mu1_hat+sqrt((p*(n-1)/(n-p))*finv(1-alpha,p,n-p))*sqrt(Sigma_hat(1,1)/n)]
line([mu1_sim_CI(1) mu1_sim_CI(1)],[650 1000],'LineWidth',2,'Color','r')
line([mu1_sim_CI(2) mu1_sim_CI(2)],[650 1000],'LineWidth',2,'Color','r')
mu2_sim_CI = [mu2_hat-sqrt((p*(n-1)/(n-p))*finv(1-alpha,p,n-p))*sqrt(Sigma_hat(2,2)/n) mu2_hat+sqrt((p*(n-1)/(n-p))*finv(1-alpha,p,n-p))*sqrt(Sigma_hat(2,2)/n)]
line([1650 2050],[mu2_sim_CI(1) mu2_sim_CI(1)],'LineWidth',2,'Color','r')
line([1650 2050],[mu2_sim_CI(2) mu2_sim_CI(2)],'LineWidth',2,'Color','r')
%--------------------------------------------------------------------------
% Marginal confidence intervals for mu (not for observations !)
% This is narrow / optimistic
%--------------------------------------------------------------------------
mu1_marg_CI = [mu1_hat-tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(1,1)/n) mu1_hat+tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(1,1)/n)]
line([mu1_marg_CI(1) mu1_marg_CI(1)],[650 1000],'LineWidth',2,'Color','b')
line([mu1_marg_CI(2) mu1_marg_CI(2)],[650 1000],'LineWidth',2,'Color','b')
mu2_marg_CI = [mu2_hat-tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(2,2)/n) mu2_hat+tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(2,2)/n)]
line([1650 2050],[mu2_marg_CI(1) mu2_marg_CI(1)],'LineWidth',2,'Color','b')
line([1650 2050],[mu2_marg_CI(2) mu2_marg_CI(2)],'LineWidth',2,'Color','b')
%--------------------------------------------------------------------------
% Bonferroni simultaneous confidence intervals for mu (not for observations !)
% This is just right :)
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
% Exact test of hypothesis H0: mu = mu0
%--------------------------------------------------------------------------
figure(6)
axis off
% % My hypothesis is that mu_1 ~= 1750 and mu_2 ~= 950
mu0 = [1750 950];
text(0,0.8,['Exact test of hypothesis   H_0: {\mu} = {\mu}_0 = ' mat2str(mu0) ',  {\alpha} = ' num2str(alpha)],'Fontsize',16)
% % This variable calculates the Hotelling's T^2 statistic.
% % It measures how far the sample means are from the specified means,
% % normalized by the covariance structure of the data. Essentially,
% % this normalization helps to ensure that the test is sensitive to the structure of the data.
T2 = (mu_hat-mu0)*inv(Sigma_hat/n)*(mu_hat-mu0)'
% % The critical value is a threshold derived from the F-distribution that
% % determines the cutoff point for deciding whether to reject the null
% % hypothesis. It is calculated based on the significance level α (e.g., 0.05)
% % and the degrees of freedom of the test, which depend on the number of
% % variables p and the sample size n. It's called a null hypothesis because
% % that indicates a LACK of relationship between the data.
critical_value = (p*(n-1)/(n-p))*finv(1-alpha,p,n-p)
text(0,0.6,['T^2 =  ' num2str(T2) '  >  critical value  = p(n-1)/(n-p)*F_p_,_n_-_p({\alpha})  =  ' num2str(critical_value)],'Fontsize',14)
% % This is the alternate way of doing it. The p-value is the probability
% % of observing a test statistic as extreme as, or more extreme than, the
% % one calculated from your sample data, assuming that the null hypothesis
% % is true. It quantifies the strength of evidence against the null hypothesis.
p_value = 1-fcdf(((n-p)/(p*(n-1)))*T2,p,n-p);
text(0,0.45,['p-value =  ' num2str(p_value)],'Fontsize',14)