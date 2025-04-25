clc;clear;close all
format compact
%--------------------------------------------------------------------------
disp('-----------------------------------------------------')
SIGMA = [4 1 2; 1 9 -3; 2 -3 25]
[P,LAMBDA] = eig(SIGMA)
disp('-----------------------------------------------------')
P_LAMBDA_PT = P*LAMBDA*P'
disp('-----------------------------------------------------')
detP = det(P)
P_PT = P*P'
disp('-----------------------------------------------------')
detSIGMA = det(SIGMA)
eigenvalues_product = prod(diag(LAMBDA))
disp('-----------------------------------------------------')
traceSIGMA = trace(SIGMA)
eigenvalues_sum = sum(diag(LAMBDA))
disp('-----------------------------------------------------')
sqrtSIGMA = P*diag(sqrt(diag(LAMBDA)))*P'
sqr_sqrtSIGMA = sqrtSIGMA^2
disp('-----------------------------------------------------')
reciprocal_sqrtSIGMA = P*diag(1./sqrt(diag(LAMBDA)))*P'
sqrtSIGMA_reciprocal_sqrtSIGMA = sqrtSIGMA*reciprocal_sqrtSIGMA
disp('-----------------------------------------------------')
command_sqrtSIGMA = sqrt(SIGMA)
command_sqrtmSIGMA = sqrtm(SIGMA)
disp('-----------------------------------------------------')
%--------------------------------------------------------------------------
