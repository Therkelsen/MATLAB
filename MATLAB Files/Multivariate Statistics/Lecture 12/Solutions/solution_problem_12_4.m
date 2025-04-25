clear;close all;clc;
format compact

%----------------------------------------------------
% load data
%----------------------------------------------------
load('M:\Statistik\MULTISTAT\Lektion 12  Cluster Analysis\Course material\dataset_problem_12_4')
disp('---------------------------------------------------------------------')
[n, p] = size(X)
disp('---------------------------------------------------------------------')
load('M:\Statistik\MULTISTAT\Lektion 12  Cluster Analysis\Course material\city_names_problem_12_4')

%----------------------------------------------------
% K = 2
%----------------------------------------------------
figure
subplot(7,8,1:8)
axis off
text(0,1.2,'Clustering with Finite Gaussian Mixture Model, example','Fontsize',28)
text(0,0.55,'GMM model for observations, X:   f_X(x) = {\Sigma}_k_=_1_._._K[p_kf_Y_k(x)],  x = [x_1 ... x_p]^T, where  Y_k ~ N_p({\mu}_k,{\Sigma}_k),  {\Sigma}_k_=_1_._._K[p_k] = 1,  p_k {\geq} 0','Fontsize',18)
subplot(7,8,11:8:51)
axis off
line([0.5 0.5],[0 1],'Color','k','LineWidth',3)
subplot(7,8,14:8:54)
axis off
line([0.5 0.5],[0 1],'Color','k','LineWidth',3)
subplot(7,8,[17 18 25 26])
plot(X(:,1),X(:,2),'k*','LineWidth',1),xlim([0 20]),ylim([0 60]),grid,axis square,xlabel('x_1','Fontsize',14),ylabel('x_2','Fontsize',14)
subplot(7,8,9:10)
axis off
text(0.23,0.8,'Observations','Fontsize',20)
text(0.3,-0.1,['n = ' num2str(n) , ', p = ' num2str(p)],'Fontsize',14)
%----------------------------------------------------
% K = 2: estimate parameters for MVN's
%----------------------------------------------------
disp('  ')
disp('---------------------------------------------------------------------')
K = 2
disp('---------------------------------------------------------------------')
obj2 = gmdistribution.fit(X,K,'Start','randSample','Replicates',100);
estimated_priors = obj2.PComponents
disp('---------------------------------------------------------------------')
estimated_mean_vectors = obj2.mu
disp('---------------------------------------------------------------------')
estimated_covariance_matrices = obj2.Sigma
disp('---------------------------------------------------------------------')
%----------------------------------------------------
% K = 2: plot estimated pdf and contour
%----------------------------------------------------
p1_hat = obj2.PComponents(1);
mu1_hat = obj2.mu(1,:);
Sigma1_hat = obj2.Sigma(:,:,1);
%---------------------------
p2_hat = obj2.PComponents(2);
mu2_hat = obj2.mu(2,:);
Sigma2_hat = obj2.Sigma(:,:,2);
%----------------------------
x1 = 0:0.5:20;
x2 = 0:1.5:60;
[X1,X2] = meshgrid(x1,x2);
F2 = p1_hat*mvnpdf([X1(:) X2(:)],mu1_hat,Sigma1_hat) + p2_hat*mvnpdf([X1(:) X2(:)],mu2_hat,Sigma2_hat);
F2 = reshape(F2,length(x2),length(x1));
%------------------------------
subplot(7,8,[20 21 28 29])
surf(x1,x2,F2);
caxis([min(F2(:))-.5*range(F2(:)),max(F2(:))]);
xlabel('x_1','Fontsize',14); ylabel('x_2','Fontsize',14)
%--------------------------------
subplot(7,8,[44 45 52 53])
contour(x1,x2,F2),axis square,axis([0 20 0 60]),xlabel('x_1','Fontsize',14); ylabel('x_2','Fontsize',14)
subplot(7,8,12:13)
axis off
text(0,0.8,['K = ' num2str(K) ':  Estimation of GMM'],'Fontsize',20)
text(-0.15,0.45,'(see estimated parameters p_k, {\mu}_k, {\Sigma}_k in Command Window)','Fontsize',12)
text(0.2,-0.1,'estimated GMM pdf','Fontsize',14)
subplot(7,8,36:37)
axis off
text(0,-0.1,'contour plot of estimated GMM pdf','Fontsize',14)
%----------------------------------------------------
% K = 2: assign observations to clusters using MAP
%----------------------------------------------------
idx2 = cluster(obj2,X);
idx2_cluster_1 = (find(idx2 == 1))'
n2_cluster_1 = length(idx2_cluster_1);
idx2_cluster_1_names = names(idx2_cluster_1(1:end))
disp('---------------------------------------------------------------------')
idx2_cluster_2 = (find(idx2 == 2))'
n2_cluster_2 = length(idx2_cluster_2);
idx2_cluster_2_names = names(idx2_cluster_2(1:end))
disp('---------------------------------------------------------------------')
K2_cluster_1 = X(idx2_cluster_1,:);
K2_cluster_2 = X(idx2_cluster_2,:);
subplot(7,8,[23 24 31 32])
plot(K2_cluster_1(:,1),K2_cluster_1(:,2),'r*','LineWidth',1),xlim([0 20]),ylim([0 60]),grid,axis square,xlabel('x_1','Fontsize',14),ylabel('x_2','Fontsize',14)
hold on
plot(K2_cluster_2(:,1),K2_cluster_2(:,2),'b*','LineWidth',1)
hold off
%----------------------------------------------------
% K = 2: illustrate posteriors
%----------------------------------------------------
posterior_probs = posterior(obj2,X);
K2_posteriors_cluster_1 = posterior_probs(idx2_cluster_1,:);
K2_posteriors_cluster_2 = posterior_probs(idx2_cluster_2,:);
subplot(7,8,[47 48 55 56])
plot(1:n2_cluster_1,K2_posteriors_cluster_1(:,1),'r','LineWidth',3),xlim([0 n+1]),ylim([-0.1 1.1]),grid,axis square,xlabel('observation','Fontsize',14),ylabel('probability','Fontsize',14)
hold on
plot(1:n2_cluster_1,K2_posteriors_cluster_1(:,2),'b','LineWidth',1)
plot(n2_cluster_1+1:n,K2_posteriors_cluster_2(:,1),'r','LineWidth',1)
plot(n2_cluster_1+1:n,K2_posteriors_cluster_2(:,2),'b','LineWidth',3)
hold off
subplot(7,8,15:16)
axis off
text(0.03,0.8,'Clustering with GMM','Fontsize',20)
text(-0.05,0.15,'assignment of observations to clusters','Fontsize',14)
text(-0.05,-0.1,'using MAP (Maximum Aposteriori Prob.)','Fontsize',14)
subplot(7,8,39:40)
axis off
text(-0.1,0.15,'posterior probabilites, P(k|x_j), k=1:K, j=1:n,','Fontsize',14)
text(-0.1,-0.1,'for each observation belonging to each cluster','Fontsize',14)

%-------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------

%----------------------------------------------------
% K = 3
%----------------------------------------------------
figure
subplot(7,8,1:8)
axis off
text(0,1.2,'Clustering with Finite Gaussian Mixture Model, example','Fontsize',28)
text(0,0.55,'GMM model for observations, X:   f_X(x) = {\Sigma}_k_=_1_._._K[p_kf_Y_k(x)],  x = [x_1 ... x_p]^T, where  Y_k ~ N_p({\mu}_k,{\Sigma}_k),  {\Sigma}_k_=_1_._._K[p_k] = 1,  p_k {\geq} 0','Fontsize',18)
subplot(7,8,11:8:51)
axis off
line([0.5 0.5],[0 1],'Color','k','LineWidth',3)
subplot(7,8,14:8:54)
axis off
line([0.5 0.5],[0 1],'Color','k','LineWidth',3)
subplot(7,8,[17 18 25 26])
plot(X(:,1),X(:,2),'k*','LineWidth',1),xlim([0 20]),ylim([0 60]),grid,axis square,xlabel('x_1','Fontsize',14),ylabel('x_2','Fontsize',14)
subplot(7,8,9:10)
axis off
text(0.23,0.8,'Observations','Fontsize',20)
text(0.3,-0.1,['n = ' num2str(n) , ', p = ' num2str(p)],'Fontsize',14)
%----------------------------------------------------
% K = 3: estimate parameters for MVN's
%----------------------------------------------------
disp('  ')
disp('---------------------------------------------------------------------')
K = 3
disp('---------------------------------------------------------------------')
obj3 = gmdistribution.fit(X,K,'Start','randSample','Replicates',100);
estimated_priors = obj3.PComponents
disp('---------------------------------------------------------------------')
estimated_mean_vectors = obj3.mu
disp('---------------------------------------------------------------------')
estimated_covariance_matrices = obj3.Sigma
disp('---------------------------------------------------------------------')
%----------------------------------------------------
% K = 3: plot estimated pdf and contour
%----------------------------------------------------
p1_hat = obj3.PComponents(1);
mu1_hat = obj3.mu(1,:);
Sigma1_hat = obj3.Sigma(:,:,1);
%---------------------------
p2_hat = obj3.PComponents(2);
mu2_hat = obj3.mu(2,:);
Sigma2_hat = obj3.Sigma(:,:,2);
%---------------------------
p3_hat = obj3.PComponents(3);
mu3_hat = obj3.mu(3,:);
Sigma3_hat = obj3.Sigma(:,:,3);
%----------------------------
x1 = 0:0.5:20;
x2 = 0:1.5:60;
[X1,X2] = meshgrid(x1,x2);
F3 = p1_hat*mvnpdf([X1(:) X2(:)],mu1_hat,Sigma1_hat) + p2_hat*mvnpdf([X1(:) X2(:)],mu2_hat,Sigma2_hat) + p3_hat*mvnpdf([X1(:) X2(:)],mu3_hat,Sigma3_hat);
F3 = reshape(F3,length(x2),length(x1));
%------------------------------
subplot(7,8,[20 21 28 29])
surf(x1,x2,F3);
caxis([min(F3(:))-.5*range(F3(:)),max(F3(:))]);
xlabel('x_1','Fontsize',14); ylabel('x_2','Fontsize',14)
%--------------------------------
subplot(7,8,[44 45 52 53])
contour(x1,x2,F3),axis square,axis([0 20 0 60]),xlabel('x_1','Fontsize',14); ylabel('x_2','Fontsize',14)
subplot(7,8,12:13)
axis off
text(0,0.8,['K = ' num2str(K) ':  Estimation of GMM'],'Fontsize',20)
text(-0.15,0.45,'(see estimated parameters p_k, {\mu}_k, {\Sigma}_k in Command Window)','Fontsize',12)
text(0.2,-0.1,'estimated GMM pdf','Fontsize',14)
subplot(7,8,36:37)
axis off
text(0,-0.1,'contour plot of estimated GMM pdf','Fontsize',14)
%----------------------------------------------------
% K = 3: assign observations to clusters using MAP
%----------------------------------------------------
idx3 = cluster(obj3,X);
idx3_cluster_1 = (find(idx3 == 1))'
n3_cluster_1 = length(idx3_cluster_1);
idx3_cluster_1_names = names(idx3_cluster_1(1:end))
disp('---------------------------------------------------------------------')
idx3_cluster_2 = (find(idx3 == 2))'
n3_cluster_2 = length(idx3_cluster_2);
idx3_cluster_1_names = names(idx3_cluster_2(1:end))
disp('---------------------------------------------------------------------')
idx3_cluster_3 = (find(idx3 == 3))'
n3_cluster_3 = length(idx3_cluster_3);
idx3_cluster_1_names = names(idx3_cluster_3(1:end))
disp('---------------------------------------------------------------------')
K3_cluster_1 = X(idx3_cluster_1,:);
K3_cluster_2 = X(idx3_cluster_2,:);
K3_cluster_3 = X(idx3_cluster_3,:);
subplot(7,8,[23 24 31 32])
plot(K3_cluster_1(:,1),K3_cluster_1(:,2),'r*','LineWidth',1),xlim([0 20]),ylim([0 60]),grid,axis square,xlabel('x_1','Fontsize',14),ylabel('x_2','Fontsize',14)
hold on
plot(K3_cluster_2(:,1),K3_cluster_2(:,2),'b*','LineWidth',1)
plot(K3_cluster_3(:,1),K3_cluster_3(:,2),'g*','LineWidth',1)
hold off
%----------------------------------------------------
% K = 3: illustrate posteriors
%----------------------------------------------------
posterior_probs = posterior(obj3,X);
K3_posteriors_cluster_1 = posterior_probs(idx3_cluster_1,:);
K3_posteriors_cluster_2 = posterior_probs(idx3_cluster_2,:);
K3_posteriors_cluster_3 = posterior_probs(idx3_cluster_3,:);
subplot(7,8,[47 48 55 56])
plot(1:n3_cluster_1,K3_posteriors_cluster_1(:,1),'r','LineWidth',3),xlim([0 n+1]),ylim([-0.1 1.1]),grid,axis square,xlabel('observation','Fontsize',14),ylabel('probability','Fontsize',14)
hold on
plot(1:n3_cluster_1,K3_posteriors_cluster_1(:,2),'b','LineWidth',1)
plot(1:n3_cluster_1,K3_posteriors_cluster_1(:,3),'g','LineWidth',1)
plot(n3_cluster_1+1:n3_cluster_1+n3_cluster_2,K3_posteriors_cluster_2(:,1),'r','LineWidth',1)
plot(n3_cluster_1+1:n3_cluster_1+n3_cluster_2,K3_posteriors_cluster_2(:,2),'b','LineWidth',3)
plot(n3_cluster_1+1:n3_cluster_1+n3_cluster_2,K3_posteriors_cluster_2(:,3),'g','LineWidth',1)
plot(n3_cluster_1+n3_cluster_2+1:n,K3_posteriors_cluster_3(:,1),'r','LineWidth',1)
plot(n3_cluster_1+n3_cluster_2+1:n,K3_posteriors_cluster_3(:,2),'b','LineWidth',1)
plot(n3_cluster_1+n3_cluster_2+1:n,K3_posteriors_cluster_3(:,3),'g','LineWidth',3)
hold off
subplot(7,8,15:16)
axis off
text(0.03,0.8,'Clustering with GMM','Fontsize',20)
text(-0.05,0.15,'assignment of observations to clusters','Fontsize',14)
text(-0.05,-0.1,'using MAP (Maximum Aposteriori Prob.)','Fontsize',14)
subplot(7,8,39:40)
axis off
text(-0.1,0.15,'posterior probabilites, P(k|x_j), k=1:K, j=1:n,','Fontsize',14)
text(-0.1,-0.1,'for each observation belonging to each cluster','Fontsize',14)

%-------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------

%----------------------------------------------------
% Comparison
%----------------------------------------------------
figure
subplot(7,8,1:8)
axis off
text(-0.05,1.2,'Comparison of GMM clustering for different numbers of clusters, K, using','Fontsize',28)
text(-0.05,0.6,'Akaike Information Criterion,     AIC = 2*log L_m_a_x - 2*N_p_a_r(K)','Fontsize',18)
text(-0.05,0.1,'Bayesian Information Criterion, BIC = 2*log L_m_a_x - N_p_a_r(K)*log n','Fontsize',18)
text(0.55,0.6,'where   L_m_a_x = max L = max {\Pi}_j_=_1_:_nf_X(x_j | p_k,{\mu}_k,{\Sigma}_k)      (maximized likelihood of observations)','Fontsize',12)
text(0.55,0.2,'and       N_p_a_r(K) = K*[1 + p(p+3)/2]-1                           (# of estimated parameters in GMM-model)','Fontsize',12)
subplot(7,8,11:8:51)
axis off
line([0.5 0.5],[0 1],'Color','k','LineWidth',3)
subplot(7,8,14:8:54)
axis off
line([0.5 0.5],[0 1],'Color','k','LineWidth',3)
subplot(7,8,9:10)
axis off
text(0.34,0.7,'K = 2','Fontsize',28)
text(0.15,0,'clustering and posteriors','Fontsize',14)
subplot(7,8,12:13)
axis off
text(0.34,0.7,'K = 3','Fontsize',28)
text(0.15,0,'clustering and posteriors','Fontsize',14)
subplot(7,8,[17 18 25 26])
plot(K2_cluster_1(:,1),K2_cluster_1(:,2),'r*','LineWidth',1),xlim([0 20]),ylim([0 60]),grid,axis square
hold on
plot(K2_cluster_2(:,1),K2_cluster_2(:,2),'b*','LineWidth',1)
hold off
subplot(7,8,[33 34 41 42])
plot(1:n2_cluster_1,K2_posteriors_cluster_1(:,1),'r','LineWidth',3),xlim([0 n+1]),ylim([-0.1 1.1]),grid,axis square,xlabel('observation','Fontsize',10),ylabel('probability','Fontsize',10)
hold on
plot(1:n2_cluster_1,K2_posteriors_cluster_1(:,2),'b','LineWidth',1)
plot(n2_cluster_1+1:n,K2_posteriors_cluster_2(:,1),'r','LineWidth',1)
plot(n2_cluster_1+1:n,K2_posteriors_cluster_2(:,2),'b','LineWidth',3)
hold off
subplot(7,8,[20 21 28 29])
plot(K3_cluster_1(:,1),K3_cluster_1(:,2),'r*','LineWidth',1),xlim([0 20]),ylim([0 60]),grid,axis square
hold on
plot(K3_cluster_2(:,1),K3_cluster_2(:,2),'b*','LineWidth',1)
plot(K3_cluster_3(:,1),K3_cluster_3(:,2),'g*','LineWidth',1)
hold off
subplot(7,8,[36 37 44 45])
plot(1:n3_cluster_1,K3_posteriors_cluster_1(:,1),'r','LineWidth',3),xlim([0 n+1]),ylim([-0.1 1.1]),grid,axis square,xlabel('observation','Fontsize',10),ylabel('probability','Fontsize',10)
hold on
plot(1:n3_cluster_1,K3_posteriors_cluster_1(:,2),'b','LineWidth',1)
plot(1:n3_cluster_1,K3_posteriors_cluster_1(:,3),'g','LineWidth',1)
plot(n3_cluster_1+1:n3_cluster_1+n3_cluster_2,K3_posteriors_cluster_2(:,1),'r','LineWidth',1)
plot(n3_cluster_1+1:n3_cluster_1+n3_cluster_2,K3_posteriors_cluster_2(:,2),'b','LineWidth',3)
plot(n3_cluster_1+1:n3_cluster_1+n3_cluster_2,K3_posteriors_cluster_2(:,3),'g','LineWidth',1)
plot(n3_cluster_1+n3_cluster_2+1:n,K3_posteriors_cluster_3(:,1),'r','LineWidth',1)
plot(n3_cluster_1+n3_cluster_2+1:n,K3_posteriors_cluster_3(:,2),'b','LineWidth',1)
plot(n3_cluster_1+n3_cluster_2+1:n,K3_posteriors_cluster_3(:,3),'g','LineWidth',3)
hold off

subplot(7,8,49:50)
axis off
K = 2;
logL_2_clusters = -obj2.NlogL;
Npar_2_clusters = K*(1+p*(p+3)/2)-1;
AIC_2_clusters = -obj2.AIC;
BIC_2_clusters = -obj2.BIC;
text(0.21,0.45,['2*log L_m_a_x = ' num2str(2*logL_2_clusters,4)],'Fontsize',14)
text(0.21,0.08,['N_p_a_r = ' num2str(Npar_2_clusters)],'Fontsize',14)
text(0.21,-0.31,['AIC = ' num2str(AIC_2_clusters,4)],'Fontsize',14)
text(0.21,-0.65,['BIC = ' num2str(BIC_2_clusters,4)],'Fontsize',14)
subplot(7,8,52:53)
axis off
K = 3;
logL_3_clusters = -obj3.NlogL;
Npar_3_clusters = K*(1+p*(p+3)/2)-1;
AIC_3_clusters = -obj3.AIC;
BIC_3_clusters = -obj3.BIC;
text(0.21,0.45,['2*log L_m_a_x = ' num2str(2*logL_3_clusters,4)],'Fontsize',14)
text(0.21,0.08,['N_p_a_r = ' num2str(Npar_3_clusters)],'Fontsize',14)
text(0.21,-0.31,['AIC = ' num2str(AIC_3_clusters,4) '  (K=3 best AIC)'],'Fontsize',14)
text(0.21,-0.65,['BIC = ' num2str(BIC_3_clusters,4) '     (K=3 best BIC)'],'Fontsize',14)

%---------------------------------------------------------------------------------------------------------------------------------
