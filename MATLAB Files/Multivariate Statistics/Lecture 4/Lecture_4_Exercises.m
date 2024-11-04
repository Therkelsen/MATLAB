clc; clear; close all;
format compact

utils = Utils;

disp("Problem 4.1")

%--------------------------------------------------------------------------
% load data
%--------------------------------------------------------------------------
data = load("Lecture 4/dataset_problem_4_1.mat")
disp('-------------------------------------------------------------------')
disp('Bartlett test,  H0: SIGMA_female = SIGMA_male')
disp('-------------------------------------------------------------------')
%--------------------------------------------------------------------------
% descriptive statistics
%--------------------------------------------------------------------------
% % Calculate descriptive statustics (Mean vector, Covariance matrix, range vector)
[mu_hat_male, Sigma_hat_male, range_male] = utils.calculate_descriptive_statistics(data.X_male, true, true, "Male turtles");
% % Column 1 and 2 of the mean vector
mu1_hat_male = mu_hat_male(1); mu2_hat_male = mu_hat_male(2);
% % Standard deviations of diagonals of covariance matrix (std = sqrt(cov))
sigma1_hat_male = sqrt(Sigma_hat_male(1,1)); sigma2_hat_male = sqrt(Sigma_hat_male(2,2));
% % The biggest covariance element of the diagonals
maxvar_male = max(Sigma_hat_male(1,1),Sigma_hat_male(2,2));

% % Plot the observations using a scatter plot
% % n rows p columns, each p is a variable
[n_male, p_male] = size(data.X_male)

% % Calculate descriptive statustics (Mean vector, Covariance matrix, range vector)
[mu_hat_female, Sigma_hat_female, range_female] = utils.calculate_descriptive_statistics(data.X_female, true, true, "Female turtles");
% % Column 1 and 2 of the mean vector
mu1_hat_female = mu_hat_female(1); mu2_hat_female = mu_hat_female(2);
% % Standard deviations of diagonals of covariance matrix (std = sqrt(cov))
sigma1_hat_female = sqrt(Sigma_hat_female(1,1)); sigma2_hat_female = sqrt(Sigma_hat_female(2,2));
% % The biggest covariance element of the diagonals
maxvar_female = max(Sigma_hat_female(1,1),Sigma_hat_female(2,2));

% % Plot the observations using a scatter plot
% % n rows p columns, each p is a variable
[n_female, p_female] = size(data.X_female)

if (p_male == p_female)
    p = p_male
end

Sigma_hat_pooled = ((n_female-1)*Sigma_hat_female + (n_male-1)*Sigma_hat_male)/(n_female+n_male-2)

% We wish to test if the covariance matrices for female and male
% turtles can be said to be equal, that is to test
% H0: Cov_female = Cov_male against H1: Cov_female != Cov_male

%--------------------------------------------------------------------------
%  Bartlett test (also known as Box's M-test) for equal covariance matrices
%  This test is used to determine if two groups (female and male) have 
%  statistically similar covariance matrices. If the null hypothesis is rejected,
%  we conclude that the covariance matrices are significantly different.
%--------------------------------------------------------------------------

% Calculate the main part of the test statistic (T) using the log determinants 
% of the pooled and group-specific covariance matrices. 
% 'Sigma_hat_pooled' is the pooled covariance matrix.
T = (n_female + n_male - 2) * log(det(Sigma_hat_pooled)) ...
    - (n_female - 1) * log(det(Sigma_hat_female)) ...
    - (n_male - 1) * log(det(Sigma_hat_male));

% Apply a correction factor to adjust the test statistic for sample sizes.
% This correction accounts for small sample bias in multivariate data and 
% depends on the number of variables (p) and sample sizes of each group.
correction_factor = 1 - ((2 * p^2 + 3 * p - 1) / (6 * (p + 1))) ...
                    * (1 / (n_female - 1) + 1 / (n_male - 1) - 1 / (n_female + n_male - 2));

% Calculate the corrected test statistic by applying the correction factor to T.
test_statistic = correction_factor * T

% Define significance level (alpha) for the hypothesis test
alpha = 0.05;

% Calculate the degrees of freedom for the chi-square distribution.
% Degrees of freedom depend on the number of variables in the covariance matrices.
df = p * (p + 1) / 2;

% Determine the critical value from the chi-square distribution.
% If the test statistic exceeds this value, we reject the null hypothesis.
critical_value = chi2inv(1 - alpha, df)

% Evaluate whether to reject the null hypothesis (H0) of equal covariance matrices.
% If 'reject_H0' is true, we conclude that the groups have significantly different covariances.
reject_H0 = test_statistic > critical_value

% Calculate the p-value to determine the exact significance level of the test statistic.
p_value = 1 - chi2cdf(test_statistic, df)

% Display a separation line in the output for readability.
disp('-------------------------------------------------------------------')
%------------------------------------------------------------------------------------------
% Model check for MALE data
%------------------------------------------------------------------------
mahalanobis_distances_male = utils.calculate_mahalanobis_distances(data.X_male, mu_hat_male, Sigma_hat_male, true, true, "Male turtles");
%------------------------------------------------------------------------------------------
% Model check for FEMALE data
%------------------------------------------------------------------------
mahalanobis_distances_female = utils.calculate_mahalanobis_distances(data.X_female, mu_hat_female, Sigma_hat_female, true, true, "Female turtles");

%%
clc; clear; close all;
format compact

utils = Utils;

disp("Problem 4.2")

%--------------------------------------------------------------------------
% load data
%--------------------------------------------------------------------------
data = load("Lecture 4/dataset_problem_4_2.mat")
disp('-------------------------------------------------------------------')
disp('Bartlett test,  H0: SIGMA_female = SIGMA_male')
disp('-------------------------------------------------------------------')
%--------------------------------------------------------------------------
% descriptive statistics
%--------------------------------------------------------------------------
% % Calculate descriptive statustics (Mean vector, Covariance matrix, range vector)
[mu_hat_male, Sigma_hat_male, range_male] = utils.calculate_descriptive_statistics(data.X_male, true, true, "Male birds");
% % Column 1 and 2 of the mean vector
mu1_hat_male = mu_hat_male(1); mu2_hat_male = mu_hat_male(2);
% % Standard deviations of diagonals of covariance matrix (std = sqrt(cov))
sigma1_hat_male = sqrt(Sigma_hat_male(1,1)); sigma2_hat_male = sqrt(Sigma_hat_male(2,2));
% % The biggest covariance element of the diagonals
maxvar_male = max(Sigma_hat_male(1,1),Sigma_hat_male(2,2));

% % Plot the observations using a scatter plot
% % n rows p columns, each p is a variable
[n_male, p_male] = size(data.X_male)

% % Calculate descriptive statustics (Mean vector, Covariance matrix, range vector)
[mu_hat_female, Sigma_hat_female, range_female] = utils.calculate_descriptive_statistics(data.X_female, true, true, "Female birds");
% % Column 1 and 2 of the mean vector
mu1_hat_female = mu_hat_female(1); mu2_hat_female = mu_hat_female(2);
% % Standard deviations of diagonals of covariance matrix (std = sqrt(cov))
sigma1_hat_female = sqrt(Sigma_hat_female(1,1)); sigma2_hat_female = sqrt(Sigma_hat_female(2,2));
% % The biggest covariance element of the diagonals
maxvar_female = max(Sigma_hat_female(1,1),Sigma_hat_female(2,2));

% % Plot the observations using a scatter plot
% % n rows p columns, each p is a variable
[n_female, p_female] = size(data.X_female)

if (p_male == p_female)
    p = p_male
end

if (n_male == n_female)
    n = n_male
end

Sigma_hat_pooled = ((n_female-1)*Sigma_hat_female + (n_male-1)*Sigma_hat_male)/(n_female+n_male-2)

% We wish to test if the covariance matrices for female and male
% birds can be said to be equal, that is to test
% H0: Cov_female = Cov_male against H1: Cov_female != Cov_male

%--------------------------------------------------------------------------
%  Bartlett test (also known as Box's M-test) for equal covariance matrices
%  This test is used to determine if two groups (female and male) have 
%  statistically similar covariance matrices. If the null hypothesis is rejected,
%  we conclude that the covariance matrices are significantly different.
%--------------------------------------------------------------------------

% Calculate the main part of the test statistic (T) using the log determinants 
% of the pooled and group-specific covariance matrices. 
% 'Sigma_hat_pooled' is the pooled covariance matrix.
T = (n_female + n_male - 2) * log(det(Sigma_hat_pooled)) ...
    - (n_female - 1) * log(det(Sigma_hat_female)) ...
    - (n_male - 1) * log(det(Sigma_hat_male));

% Apply a correction factor to adjust the test statistic for sample sizes.
% This correction accounts for small sample bias in multivariate data and 
% depends on the number of variables (p) and sample sizes of each group.
correction_factor = 1 - ((2 * p^2 + 3 * p - 1) / (6 * (p + 1))) ...
                    * (1 / (n_female - 1) + 1 / (n_male - 1) - 1 / (n_female + n_male - 2));

% Calculate the corrected test statistic by applying the correction factor to T.
test_statistic = correction_factor * T

% Define significance level (alpha) for the hypothesis test
alpha = 0.05;

% Calculate the degrees of freedom for the chi-square distribution.
% Degrees of freedom depend on the number of variables in the covariance matrices.
df = p * (p + 1) / 2;

% Determine the critical value from the chi-square distribution.
% If the test statistic exceeds this value, we reject the null hypothesis.
critical_value = chi2inv(1 - alpha, df)

% Evaluate whether to reject the null hypothesis (H0) of equal covariance matrices.
% If 'reject_H0' is true, we conclude that the groups have significantly different covariances.
reject_H0 = test_statistic > critical_value

% Calculate the p-value to determine the exact significance level of the test statistic.
p_value = 1 - chi2cdf(test_statistic, df)

% Display a separation line in the output for readability.
disp('-------------------------------------------------------------------')
%------------------------------------------------------------------------------------------
% Model check for MALE data
%------------------------------------------------------------------------
mahalanobis_distances_male = utils.calculate_mahalanobis_distances(data.X_male, mu_hat_male, Sigma_hat_male, true, true, "Male birds");
%------------------------------------------------------------------------------------------
% Model check for FEMALE data
%------------------------------------------------------------------------
mahalanobis_distances_female = utils.calculate_mahalanobis_distances(data.X_female, mu_hat_female, Sigma_hat_female, true, true, "Female birds");

% We wish to test for equality of mean vectors for the male
% and female populations of hook-billed kites, that is to test
% H0mu: Mu_female = Mu_male against H1mu: Mu_female != Mu_male

%--------------------------------------------------------------------------
% Exact test of hypothesis H0: mu_female = mu_male
%--------------------------------------------------------------------------
% Calculate Hotelling's T-squared statistic (T2) for testing the equality 
% of mean vectors. This statistic considers the difference between the 
% sample mean vectors (mu_hat_female and mu_hat_male) and scales it by the 
% pooled covariance matrix (Sigma_hat_pooled) adjusted for the sample sizes.
T2 = (mu_hat_female - mu_hat_male) * inv(Sigma_hat_pooled * (1 / n_female + 1 / n_male)) ...
     * (mu_hat_female - mu_hat_male)'

% Compute the critical value based on the F-distribution.
% This value is used as a threshold for determining if the test statistic 
% T2 is significant at the specified alpha level (e.g., 0.05).
critical_value = (p * (n_female + n_male - 2) / (n_female + n_male - p - 1)) ...
                 * finv(1 - alpha, p, n_female + n_male - p - 1)

% Calculate the p-value for the T-squared test.
% This p-value tells us the probability of observing a test statistic as 
% extreme as T2 under the null hypothesis.
p_value = 1 - fcdf(((n_female + n_male - p - 1) / (p * (n_female + n_male - 2))) * T2, ...
                   p, n_female + n_male - p - 1)

reject_H0_exact_critical = T2 > critical_value;

%--------------------------------------------------------------------------
% Large sample approximate test of hypothesis H0: mu_female = mu_male
%--------------------------------------------------------------------------
% For large samples, we can approximate the test using the chi-square 
% distribution instead of the F-distribution. This is generally used 
% when the sample sizes are large enough to rely on asymptotic properties.
critical_value = chi2inv(1 - alpha, p)

% Calculate the approximate p-value for the large sample test.
% The p-value is based on the chi-square distribution with 'p' degrees of freedom.
p_value = 1 - chi2cdf(T2, p)

reject_H0_large_critical = T2 > critical_value;

reject_H0_large_pval = p_value < alpha;

reject_H0 = reject_H0_exact_critical && reject_H0_large_critical && reject_H0_large_pval

%--------------------------------------------------------------------------
% Further descriptive statistics for plotting the confidence region
%--------------------------------------------------------------------------
% Calculate the estimated difference between the mean vectors for each 
% variable (mu_diff_hat) for plotting purposes.
mu_diff_hat = mu_hat_female - mu_hat_male

% Extract individual components of the mean difference vector for each variable.
mu_diff1_hat = mu_diff_hat(1); % Difference in the first variable
mu_diff2_hat = mu_diff_hat(2); % Difference in the second variable

% Calculate standard deviations for each variable from the pooled covariance matrix.
sigma1_hat = sqrt(Sigma_hat_pooled(1,1)); % Standard deviation of the first variable
sigma2_hat = sqrt(Sigma_hat_pooled(2,2)); % Standard deviation of the second variable

% Determine the maximum variance of the pooled covariance matrix for scaling plots.
maxvar = max(Sigma_hat_pooled(1,1), Sigma_hat_pooled(2,2));




%--------------------------------------------------------------------------
% Confidence region for mu_female - mu_male (not for observations !)
%--------------------------------------------------------------------------
figure()
y1 = mu_diff1_hat-sqrt(Sigma_hat_pooled(1,1)):sqrt(Sigma_hat_pooled(1,1))/50:mu_diff1_hat+sqrt(Sigma_hat_pooled(1,1)); 
y2 = mu_diff2_hat-sqrt(Sigma_hat_pooled(2,2)):sqrt(Sigma_hat_pooled(2,2))/50:mu_diff2_hat+sqrt(Sigma_hat_pooled(2,2));
[Y1,Y2] = meshgrid(y1,y2);
F = mvnpdf([Y1(:) Y2(:)],mu_diff_hat,Sigma_hat_pooled);
F = reshape(F,length(y2),length(y1));
alpha = 0.05;
alpha_contour = [0.05 0.05];
c = (p*(n_female+n_male-2)/(n_female+n_male-p-1))*(1/n_female+1/n_male)*finv(1-alpha_contour,p,n_female+n_male-p-1);
C = (1/(2*pi*sqrt(det(Sigma_hat_pooled))))*exp(-c/2);
plot2d_CR_for_difference_in_mu_ellipsis (mu_diff_hat', Sigma_hat_pooled, alpha, n_female, n_male)
title({'95% Confidence region and intervals for ({\mu}_f_e_m_a_l_e_,_1-{\mu}_m_a_l_e_,_1, {\mu}_f_e_m_a_l_e_,_2-{\mu}_m_a_l_e_,_2 based on observations:';'CR (black), simult.CI (red), marg.CI (blue), Bonf.CI (green)'},'Fontsize',14)
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