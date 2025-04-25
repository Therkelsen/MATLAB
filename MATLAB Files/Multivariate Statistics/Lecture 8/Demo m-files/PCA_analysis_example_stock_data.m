clear;close all;clc;
format compact

load('dataset_stock_example');

disp('---------------------------------------------------------------------------------')
disp('Descriptive statistics estimators for stock data:')
disp('(weekly rates of return during 103 weeks for 3 banks and 2 oil companies)')
disp('---------------------------------------------------------------------------------')
n = size(X,1)
p = size(X,2)
Xmean = mean(X,1)
S = cov(X)
R = corr(X)
disp('---------------------------------------------------------------------------------')
disp('Figure 1: boxplot')
figure('Name','Box plot for original stock data')
boxplot(X,'labels',{'Bank#1','Bank#2','Bank#3','OilCompany#1','OilCompany#2'}),title('Box plot for original stock data (weekly rates of return during 103 weeks)','Fontsize',14)

disp('---------------------------------------------------------------------------------')
disp('Figure 2: scatterplot')
figure('Name','Scatter plot for original stock data')
plotmatrix(X),title('Scatter plot for original stock data (weekly rates of return during 103 weeks)','Fontsize',14)
xlabel('Bank#1                               Bank#2                                     Bank#3                                OilCompany#1                            OilCompany#2','Fontsize',14)
ylabel('OilCompany#2    OilCompany#1   Bank#3          Bank#2           Bank#1   ','Fontsize',14)

%--------------------------------------------------------------------------------------------------------------------------

disp('---------------------------------------------------------------------------------')
disp('PCA on S using ''SVD(S)''')
disp('---------------------------------------------------------------------------------')
[E_S LAMBDA_S E_S] = svd(S);
E_S = sign(E_S(1,1))*E_S;
LAMBDA_S,E_S
disp('---------------------------------------------------------------------------------')
disp('PCA on S using ''PCA(X)''')
disp('---------------------------------------------------------------------------------')
[COEFFS,SCORES,VARIANCES] = pca(X);
VARIANCES,COEFFS

disp('---------------------------------------------------------------------------------')
disp('Variance explained by principal components')
disp('---------------------------------------------------------------------------------')
PC_S_percent_explained_by_PCs = 100*diag(LAMBDA_S)/sum(diag(LAMBDA_S))
disp('---------------------------------------------------------------------------------')
disp('Figure 3: PCA(S): Pareto plot and scree plot')
figure('Name','PCA(S): Pareto plot and scree plot for stock data PC''s')
bar(1:p,PC_S_percent_explained_by_PCs),xlabel('PC'),ylabel('% variance of PC''s'), hold on
title('PCA(S): Pareto plot and scree plot for stock data PC''s','Fontsize',14)
accumulated = zeros(1,p);
accumulated(1) = PC_S_percent_explained_by_PCs(1);
for i = 2:p
    accumulated(i) = sum(PC_S_percent_explained_by_PCs(1:i));
end
plot(1:p,accumulated),grid
hold off

disp('---------------------------------------------------------------------------------')
disp('Loadings for PC 1 and 2')
disp('---------------------------------------------------------------------------------')
V = diag(diag(S));
PC_1_S_loadings_X_variables = E_S(:,1)'
PC_2_S_loadings_X_variables = E_S(:,2)'

disp('---------------------------------------------------------------------------------')
disp('Figure 4: PCA(S): Scatterplot of scores (mean-corrected) for PC 1 and 2')
% PC_1_S_scores_X_variables = X*PC_1_S_loadings_X_variables';
% PC_2_S_scores_X_variables = X*PC_2_S_loadings_X_variables';
PC_1_S_scores_X_variables = (X - ones(n,1)*mean(X))*PC_1_S_loadings_X_variables';
PC_2_S_scores_X_variables = (X - ones(n,1)*mean(X))*PC_2_S_loadings_X_variables';
figure('Name','PCA(S): Scatterplot of scores (mean-corrected) for PC 1 and 2')
plot(PC_1_S_scores_X_variables,PC_2_S_scores_X_variables,'r+','Linewidth',5),grid,axis equal  
title('PCA(S): Scatterplot of scores (mean-corrected) for PC 1 and 2','Fontsize',14)
xlabel('1st Principal Component scores')
ylabel('2nd Principal Component scores')

disp('---------------------------------------------------------------------------------')
disp('Figure 5: PCA(S): Bi-plot for PC 1 and 2 scores (mean-corrected) and loadings')
figure('Name','PCA(S): Bi-plot for PC 1 and 2 scores (mean-corrected) and loadings') 
names = char('Bank#1','Bank#2','Bank#3','OilCompany#1','OilCompany#2');
biplot([PC_1_S_loadings_X_variables' PC_2_S_loadings_X_variables'],'Linewidth',3,'scores',[PC_1_S_scores_X_variables PC_2_S_scores_X_variables],'varlabels',names)
title('PCA(S): Bi-plot for PC 1 and 2 scores (mean-corrected) and loadings (note: BIPLOT sign convention inverts sign of PC2 loadings)','Fontsize',14)

disp('---------------------------------------------------------------------------------')
disp('Correlation between PC 1 and 2 and original variables')
disp('---------------------------------------------------------------------------------')
Correlations_PC_1_S_versus_X_variables = sqrt(LAMBDA_S(1,1))*PC_1_S_loadings_X_variables./(sqrt(diag(S)))'
Correlations_PC_2_S_versus_X_variables = sqrt(LAMBDA_S(2,2))*PC_2_S_loadings_X_variables./(sqrt(diag(S)))'
disp('---------------------------------------------------------------------------------')
disp('Figure 6: PCA(S): Circle diagram for correlation between PC 1 and 2 and original variables')
figure('Name','PCA(S): Circle diagram for correlation between PC 1 and 2 and original variables') 
t = 0:.01:pi;
x = cos(t);
y = sin(t);
plot(x,y,'k','LineWidth',3),axis square,title('PCA(S): Circle diagram for correlation between PC 1 and 2 and original variables','Fontsize',14),grid,xlim([-1.2 1.2]),ylim([-1.2 1.2])
xlabel('Correlation with PC 1','Fontsize',14),ylabel('Correlation with PC 2','Fontsize',14)
hold on
plot(x,-y,'k','LineWidth',3);
plot(Correlations_PC_1_S_versus_X_variables,Correlations_PC_2_S_versus_X_variables,'r+','Linewidth',10)
text(0.15,-0.8,'Bank#1','Fontsize',12)
text(0.6,-0.72,'Bank#2','Fontsize',12)
text(0.45,-0.61,'Bank#3','Fontsize',12)
text(0.41,0.24,'OilCompany#1','Fontsize',12)
text(0.41,0.31,'OilCompany#2','Fontsize',12)
text(-0.5,0.7,'Since        {\Sigma}_j_=_1_._._p r^2(X_i,PC_j) = 1','Fontsize',14)
text(-0.5,0.55,'it follows   r^2(X_i,PC_1) + r^2(X_i,PC_2) {\leq} 1','Fontsize',14)
hold off

disp('---------------------------------------------------------------------------------')
disp('Approx. simult. Bonferroni CI for variances of PC 1 and 2')
disp('---------------------------------------------------------------------------------')
alpha = 0.1
psim = 2;
alpha_bonf = alpha/psim
CI_approx_lambda_1 = LAMBDA_S(1,1)./[1+norminv(1-alpha_bonf/2)*sqrt(2/n) 1 1-norminv(1-alpha_bonf/2)*sqrt(2/n)]
CI_approx_lambda_2 = LAMBDA_S(2,2)./[1+norminv(1-alpha_bonf/2)*sqrt(2/n) 1 1-norminv(1-alpha_bonf/2)*sqrt(2/n)]
disp('---------------------------------------------------------------------------------')

disp(' ')
disp('--------------------------------------------------------------------------------------------------------------------------')
disp(' ')

%--------------------------------------------------------------------------------------------------------------------------

disp('---------------------------------------------------------------------------------')
disp('PCA on R using ''SVD(R)''')
disp('---------------------------------------------------------------------------------')
[E_R LAMBDA_R E_R] = svd(R);
E_R = sign(E_R(1,1))*E_R;
LAMBDA_R,E_R

disp('---------------------------------------------------------------------------------')
disp('Variance explained by principal components')
disp('---------------------------------------------------------------------------------')
PC_R_percent_explained_by_PCs = 100*diag(LAMBDA_R)/sum(diag(LAMBDA_R))
disp('---------------------------------------------------------------------------------')
disp('Figure 7: PCA(R): Pareto plot and scree plot')
figure('Name','PCA(R): Pareto plot and scree plot for stock data PC''s')
bar(1:p,PC_R_percent_explained_by_PCs),xlabel('PC'),ylabel('% variance of PC''s'), hold on
title('PCA(R): Pareto plot and scree plot for stock data PC''s','Fontsize',14)
accumulated = zeros(1,p);
accumulated(1) = PC_R_percent_explained_by_PCs(1);
for i = 2:p
    accumulated(i) = sum(PC_R_percent_explained_by_PCs(1:i));
end
plot(1:p,accumulated),grid
hold off

disp('---------------------------------------------------------------------------------')
disp('Loadings for PC 1 and 2')
disp('---------------------------------------------------------------------------------')
V = diag(diag(S));
PC_1_R_loadings_Z_variables = E_R(:,1)'
PC_1_R_loadings_X_variables = ((V^(-1/2))*E_R(:,1))'
PC_2_R_loadings_Z_variables = E_R(:,2)'
PC_2_R_loadings_X_variables = ((V^(-1/2))*E_R(:,2))'

disp('---------------------------------------------------------------------------------')
disp('Figure 8: PCA(R): Scatterplot of scores (mean-corrected) for PC 1 and 2')
Z = zscore(X);
PC_1_R_scores_Z_variables = Z*PC_1_R_loadings_Z_variables';
PC_1_R_scores_X_variables = (X - ones(n,1)*mean(X))*PC_1_R_loadings_X_variables';
PC_2_R_scores_Z_variables = Z*PC_2_R_loadings_Z_variables';
PC_2_R_scores_X_variables = (X - ones(n,1)*mean(X))*PC_2_R_loadings_X_variables';
figure('Name','PCA(R): Scatterplot of scores (mean-corrected) for PC 1 and 2')
plot(PC_1_R_scores_Z_variables,PC_2_R_scores_Z_variables,'r+','Linewidth',5),grid,axis equal  
title('PCA(R): Scatterplot of scores (mean-corrected) for PC 1 and 2','Fontsize',14)
xlabel('1st Principal Component scores')
ylabel('2nd Principal Component scores')

disp('---------------------------------------------------------------------------------')
disp('Figure 9: PCA(R): Bi-plot for PC 1 and 2 scores (mean-corrected) and loadings')
figure('Name','PCA(R): Bi-plot for PC 1 and 2 scores (mean-corrected) and loadings') 
names = char('Bank#1','Bank#2','Bank#3','OilCompany#1','OilCompany#2');
biplot([PC_1_R_loadings_Z_variables' PC_2_R_loadings_Z_variables'],'Linewidth',3,'scores',[PC_1_R_scores_Z_variables PC_2_R_scores_Z_variables],'varlabels',names)
title('PCA(R): Bi-plot for PC 1 and 2 scores (mean-corrected) and loadings','Fontsize',14)

disp('---------------------------------------------------------------------------------')
disp('Correlation between PC 1 and 2 and original variables')
disp('---------------------------------------------------------------------------------')
Correlations_PC_1_R_versus_Z_or_X_variables = sqrt(LAMBDA_R(1,1))*PC_1_R_loadings_Z_variables
Correlations_PC_2_R_versus_Z_or_X_variables = sqrt(LAMBDA_R(2,2))*PC_2_R_loadings_Z_variables
disp('---------------------------------------------------------------------------------')
disp('Figure 10: PCA(R): Circle diagram for correlation between PC 1 and 2 and original variables')
figure('Name','PCA(R): Circle diagram for correlation between PC 1 and 2 and original variables') 
t = 0:.01:pi;
x = cos(t);
y = sin(t);
plot(x,y,'k','LineWidth',3),axis square,title('PCA(R): Circle diagram for correlation between PC 1 and 2 and original variables','Fontsize',14),grid,xlim([-1.2 1.2]),ylim([-1.2 1.2])
xlabel('Correlation with PC 1','Fontsize',14),ylabel('Correlation with PC 2','Fontsize',14)
hold on
plot(x,-y,'k','LineWidth',3);
plot(Correlations_PC_1_R_versus_Z_or_X_variables,Correlations_PC_2_R_versus_Z_or_X_variables,'r+','Linewidth',10)
text(0.48,-0.44,'Bank#1','Fontsize',12)
text(0.58,-0.28,'Bank#2','Fontsize',12)
text(0.48,-0.37,'Bank#3','Fontsize',12)
text(0.65,0.69,'OilCompany#1','Fontsize',12)
text(0.6,0.73,'OilCompany#2','Fontsize',12)
text(-0.5,0.7,'Since        {\Sigma}_j_=_1_._._p r^2(X_i,PC_j) = 1','Fontsize',14)
text(-0.5,0.55,'it follows   r^2(X_i,PC_1) + r^2(X_i,PC_2) {\leq} 1','Fontsize',14)
hold off

disp('---------------------------------------------------------------------------------')
disp('Approx. simult. Bonferroni CI for variances of PC 1 and 2')
disp('---------------------------------------------------------------------------------')
alpha = 0.1
psim = 2;
alpha_bonf = alpha/psim
CI_approx_lambda_1 = LAMBDA_R(1,1)./[1+norminv(1-alpha_bonf/2)*sqrt(2/n) 1 1-norminv(1-alpha_bonf/2)*sqrt(2/n)]
CI_approx_lambda_2 = LAMBDA_R(2,2)./[1+norminv(1-alpha_bonf/2)*sqrt(2/n) 1 1-norminv(1-alpha_bonf/2)*sqrt(2/n)]
disp('---------------------------------------------------------------------------------')

