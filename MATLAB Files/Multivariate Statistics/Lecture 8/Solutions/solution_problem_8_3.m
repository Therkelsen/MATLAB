clear;clc;close all
format compact

disp('-------------------------------------------------------------------')
disp('1A: CALCULATION OF PCA FOR S BASED ON EVD(S)') 
disp('-------------------------------------------------------------------')
X_noncentered = rand(4,3)
[n p] = size(X_noncentered)
X = X_noncentered - ones(n,1)*mean(X_noncentered,1)
S = (1/(n-1))*X'*X
disp('-------------------------------------------------------------------')
[E_S lambda_S E_S] = svd(S);
lambda_S,E_S
% for symmetrisk kvadratisk matrix giver eig og svd samme resultat, men
% svd anvendes i stedet for eig pga sortering af egenværdier
disp('-------------------------------------------------------------------')
PCA_1_S_loading_EVD = E_S(:,1)'
PCA_1_S_scores_EVD = X*E_S(:,1)
PCA_1_S_variance_EVD = lambda_S(1,1)
PCA_1_S_variance_explained_ratio_EVD = PCA_1_S_variance_EVD/trace(lambda_S)
disp('-------------------------------------------------------------------')
%--------------------------------------------------------------------------

disp('  ')
disp('-------------------------------------------------------------------')
disp('1B: CALCULATION OF PCA FOR S BASED ON SVD(X)  (same results as 1A)') 
disp('-------------------------------------------------------------------')
Y = (1/sqrt(n-1))*X;
[U_S D_S P_S] = svd(Y,'econ');
for i = 1:p
    sign_OK = (sign(E_S(1,i)) == sign(P_S(1,i)));
    if sign_OK == 0
        P_S(:,i) = -P_S(:,i);
        U_S(:,i) = -U_S(:,i);
    end
end
U_S,D_S,P_S
disp('-------------------------------------------------------------------')
PCA_1_S_loading_SVD = P_S(:,1)'
PCA_1_S_scores_SVD = X*P_S(:,1)
PCA_1_S_variance_SVD = D_S(1,1)^2
PCA_1_S_variance_explained_ratio_SVD = PCA_1_S_variance_SVD/trace(D_S^2)
disp('-------------------------------------------------------------------')

%----------------------------------------------------------------------------------------------------

disp('  ')
disp('  ')
disp('-------------------------------------------------------------------')
disp('2A: CALCULATION OF PCA FOR R BASED ON EVD(R)') 
disp('-------------------------------------------------------------------')
V = diag(diag(S))
R = (V^(-1/2))*S*(V^(-1/2))
disp('-------------------------------------------------------------------')
[E_R lambda_R E_R] = svd(R);
lambda_R,E_R
% for symmetrisk kvadratisk matrix giver eig og svd samme resultat, men
% svd anvendes i stedet for eig pga sortering af egenværdier
disp('-------------------------------------------------------------------')
PCA_1_R_loading_EVD = ((V^(-1/2))*E_R(:,1))'
PCA_1_R_scores_EVD = X*(V^(-1/2))*E_R(:,1)
PCA_1_R_variance_EVD = lambda_R(1,1)
PCA_1_R_variance_explained_ratio_EVD = PCA_1_R_variance_EVD/trace(lambda_R)
disp('-------------------------------------------------------------------')
%--------------------------------------------------------------------------

disp('  ')
disp('-------------------------------------------------------------------')
disp('2B: CALCULATION OF PCA FOR R BASED ON SVD(Z)  (same results as 2A)') 
disp('-------------------------------------------------------------------')
Z = (1/sqrt(n-1))*X*(V^(-1/2));
[U_R D_R P_R] = svd(Z,'econ');
for i = 1:p
    sign_OK = (sign(E_R(1,i)) == sign(P_R(1,i)));
    if sign_OK == 0
        P_R(:,i) = -P_R(:,i);
        U_R(:,i) = -U_R(:,i);
    end
end
U_R,D_R,P_R
disp('-------------------------------------------------------------------')
PCA_1_R_loading_SVD = ((V^(-1/2))*P_R(:,1))'
PCA_1_R_scores_SVD = X*(V^(-1/2))*P_R(:,1)
PCA_1_R_variance_SVD = D_R(1,1)^2
PCA_1_R_variance_explained_ratio_SVD = PCA_1_R_variance_SVD/trace(D_R^2)
disp('-------------------------------------------------------------------')
%--------------------------------------------------------------------------

