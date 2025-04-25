clear;close all;clc;
format compact

disp('--------------------------------------------------------------------')
disp('PCA(SIGMA) for bivariate distribution, original X-variables')
disp('--------------------------------------------------------------------')
disp('parameters')
disp('--------------------------------------------------------------------')
sigma1 = 1
sigma2 = 10
rho = -0.4
SIGMA = [sigma1^2 rho*sigma1*sigma2; rho*sigma1*sigma2 sigma2^2]
disp('--------------------------------------------------------------------')
disp('PCA(SIGMA) using derived equations')
disp('--------------------------------------------------------------------')
lambda = (1/2)*((sigma1^2+sigma2^2) + [1 -1]*sqrt((sigma1^2-sigma2^2)^2+(2*rho*sigma1*sigma2)^2))
e = ones(2,2);
e(2,:) = (lambda-sigma1^2)/(rho*sigma1*sigma2);
e(:,1) = e(:,1)/norm(e(:,1));
e(:,2) = e(:,2)/norm(e(:,2));
e
disp('--------------------------------------------------------------------')
disp('PCA(SIGMA) using ''SVD''')
disp('--------------------------------------------------------------------')
[E LAMBDA E] = svd(SIGMA);
LAMBDA
E
disp('--------------------------------------------------------------------')
disp('Variance explained by PC 1 and PC 2')
disp('--------------------------------------------------------------------')
variance_explained_PC_SIGMA_1 = lambda(1)/sum(lambda)
variance_explained_PC_SIGMA_2 = lambda(2)/sum(lambda)
disp('--------------------------------------------------------------------')

disp(' ')
disp('--------------------------------------------------------------------')
disp('PCA(RHO) for bivariate distribution, standardized Z-variables')
disp('--------------------------------------------------------------------')
disp('parameters')
disp('--------------------------------------------------------------------')
rho
RHO = [1 rho; rho 1]
disp('--------------------------------------------------------------------')
disp('PCA(RHO) using derived equations')
disp('--------------------------------------------------------------------')
lambda = 1 + [1 -1]*rho*sign(rho)
e = (1/sqrt(2))*[1 1; [1 -1]*sign(rho)]
disp('--------------------------------------------------------------------')
disp('PCA(RHO) using ''SVD''')
disp('--------------------------------------------------------------------')
[E LAMBDA E] = svd(RHO);
LAMBDA
E
disp('--------------------------------------------------------------------')
disp('Variance explained by PC 1 and PC 2')
disp('--------------------------------------------------------------------')
variance_explained_PC_RHO_1 = (1 + rho*sign(rho))/2
variance_explained_PC_RHO_2 = (1 - rho*sign(rho))/2
disp('--------------------------------------------------------------------')