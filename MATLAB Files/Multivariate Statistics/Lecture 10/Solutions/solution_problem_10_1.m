clc;clear;close all
format compact
load('C:\Users\cv\OneDrive - Syddansk Universitet\M-drive moved by SDU-IT\Statistik\MULTISTAT\Lektion 10  Two level LDA + QDA\Course material\dataset_problem_10_1')
%------------------------------------------------------------------------
% LDA: two populations, different mu, same SIGMA
%------------------------------------------------------------------------
p = 2;
n1 = 37;
X1 = X(1:n1,2:3);
n2 = 29;
X2 = X(n1+1:end,2:3);
%------------------------------------------------------------------------
% descriptive statistics
%------------------------------------------------------------------------
x1_avg = mean(X1);
S1 = cov(X1);
x2_avg = mean(X2);
S2 = cov(X2);
%------------------------------------------------------------------------
% plot samples and estimated prediction ellipses
%------------------------------------------------------------------------
figure('Name','Two populations, different mu, same SIGMA')
subplot(4,3,1:3)
axis off
text(0,1,'Discriminant Analysis (classification) for two bivariate (p=2) populations, case I','Fontsize',28)
text(0,0.6,'Training data with known population membership:','Fontsize',20)
text(0,0.3,['Population i=1 sample: n_1 = ' num2str(n1) ' observations'],'Color','r','Fontsize',20)
text(0,0.1,['Population i=2 sample: n_2 = ' num2str(n2) ' observations'],'Color','b','Fontsize',20)
alpha = [0.5 0.1];
text(0,-0.1,['Also is shown (1-{\alpha})*100% estimated prediction ellipses for  1-{\alpha}  = ' mat2str(1-alpha) '  (assuming bivariate normal distributions)'],'Fontsize',14)
subplot(4,3,4:12)
y1 = x1_avg(1)-3*sqrt(S1(1,1)):6*sqrt(S1(1,1))/50:x1_avg(1)+3*sqrt(S1(1,1)); 
y2 = x1_avg(2)-3*sqrt(S1(2,2)):6*sqrt(S1(2,2))/50:x1_avg(2)+3*sqrt(S1(2,2));
[Y1,Y2] = meshgrid(y1,y2);
F1 = mvnpdf([Y1(:) Y2(:)],x1_avg,S1);
F1 = reshape(F1,length(y2),length(y1));
plot(X1(:,1),X1(:,2),'ro','LineWidth',3),grid
xlabel('x_i_1_j,  (i=1: j=1,...,n_1)  (i=2: j=1,...,n_2)','Fontsize',16); ylabel('x_i_2_j','Fontsize',16)
hold on
c = chi2inv(1-alpha,p);
C1 = (1/(2*pi*sqrt(det(S1)))).*exp(-c/2);
contour(y1,y2,F1,C1,'Color','r','LineWidth',1)%,axis square
%---------------------------------
z1 = x2_avg(1)-3*sqrt(S2(1,1)):6*sqrt(S2(1,1))/50:x2_avg(1)+3*sqrt(S2(1,1)); 
z2 = x2_avg(2)-3*sqrt(S2(2,2)):6*sqrt(S2(2,2))/50:x2_avg(2)+3*sqrt(S2(2,2));
[Z1,Z2] = meshgrid(z1,z2);
F2 = mvnpdf([Z1(:) Z2(:)],x2_avg,S2);
F2 = reshape(F2,length(z2),length(z1));
plot(X2(:,1),X2(:,2),'b*','LineWidth',3)
c = chi2inv(1-alpha,p);
C2 = (1/(2*pi*sqrt(det(S2)))).*exp(-c/2);
contour(z1,z2,F2,C2,'Color','b','LineWidth',1)%,axis square
hold off
%------------------------------------------------------------------------
% MVN model check for population 1 sample
%------------------------------------------------------------------------
figure('Name','MVN model check for population 1 sample')
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for observations from POPULATION 1 possibly being from bivariate normal distribution','Fontsize',20)
subplot(5,2,3)
axis off
text(0,1,'Descriptive statistics:','Fontsize',18)
text(0.5,1,['n_1 = ' num2str(n1)],'Fontsize',16)
text(0.5,0.75,['{\mu}_1_,_h_a_t = x_1 = ' mat2str(x1_avg,3)],'Fontsize',16)
text(0.635,0.96,'_','Fontsize',16),warning('off')
text(0.5,0.5,['{\Sigma}_1_,_h_a_t = S_1 = ' mat2str(S1,3)],'Fontsize',16)
subplot(5,2,5:2:9)
plotmatrix(X1,'ro')
xlabel('x_1_1_j                                          x_1_2_j','Fontsize',16)
ylabel('x_1_2_j                           x_1_1_j','Fontsize',16)
title(['Scattermatrix of observations'],'Fontsize',16)
subplot(5,2,4)
axis off
text(0,1,'Sample Mahalanobis distances, j=1..n_1 :','Fontsize',18)
text(0,0.7,'d_1_j^2 = (x_1_j - x_1)^T S_1^-^1 (x_1_j - x_1) ~ {\chi}_2^2   (approximately)','Fontsize',16)
text(0.178,0.91,'_','Fontsize',16),warning('off')
text(0.408,0.91,'_','Fontsize',16),warning('off')
subplot(5,2,6:2:10)
mu1_hat = x1_avg;
d1j_sqr = zeros(1,n1);
for j = 1:n1
    d1j_sqr(j) = (X1(j,:)-mu1_hat)*inv(S1)*(X1(j,:)-mu1_hat)';
end
df = 2;
z_j = chi2rnd(df,1,n1);
% z2 = chi2rnd(20,1,n);
qqplot(d1j_sqr,z_j),grid,xlabel('quantiles for d_1_j^2','Fontsize',16),ylabel('quantiles for {\chi}_2^2 distribution','Fontsize',16),...
    title('qq-plot for d_1_j^2 versus {\chi}_2^2 distribution','Fontsize',16)
%------------------------------------------------------------------------
% MVN model check for population 2 sample
%------------------------------------------------------------------------
figure('Name','MVN model check for population 2 sample')
subplot(5,2,1:2)
axis off
text(0,0.7,'Model check for observations from POPULATION 2 possibly being from bivariate normal distribution','Fontsize',20)
subplot(5,2,3)
axis off
text(0,1,'Descriptive statistics:','Fontsize',18)
text(0.5,1,['n_2 = ' num2str(n2)],'Fontsize',16)
text(0.5,0.75,['{\mu}_2_,_h_a_t = x_2 = ' mat2str(x2_avg,3)],'Fontsize',16)
text(0.635,0.96,'_','Fontsize',16),warning('off')
text(0.5,0.5,['{\Sigma}_2_,_h_a_t = S_2 = ' mat2str(S2,3)],'Fontsize',16)
subplot(5,2,5:2:9)
plotmatrix(X2,'b*')
xlabel('x_2_1_j                                          x_2_2_j','Fontsize',16)
ylabel('x_2_2_j                           x_2_1_j','Fontsize',16)
title(['Scattermatrix of observations'],'Fontsize',16)
subplot(5,2,4)
axis off
text(0,1,'Sample Mahalanobis distances, j=1..n_2 :','Fontsize',18)
text(0,0.7,'d_2_j^2 = (x_2_j - x_2)^T S_2^-^1 (x_2_j - x_2) ~ {\chi}_2^2   (approximately)','Fontsize',16)
text(0.178,0.91,'_','Fontsize',16),warning('off')
text(0.408,0.91,'_','Fontsize',16),warning('off')
subplot(5,2,6:2:10)
mu2_hat = x2_avg;
d2j_sqr = zeros(1,n2);
for j = 1:n2
    d2j_sqr(j) = (X2(j,:)-mu2_hat)*inv(S2)*(X2(j,:)-mu2_hat)';
end
df = 2;
z_j = chi2rnd(df,1,n2);
% z2 = chi2rnd(20,1,n);
qqplot(d2j_sqr,z_j),grid,xlabel('quantiles for d_2_j^2','Fontsize',16),ylabel('quantiles for {\chi}_2^2 distribution','Fontsize',16),...
    title('qq-plot for d_2_j^2 versus {\chi}_2^2 distribution','Fontsize',16)
%------------------------------------------------------------------------
% Bartlett test for equal covariance matrices for both populations
%------------------------------------------------------------------------
figure('Name','Bartlett and Hotelling tests --> LDA')
subplot(5,3,[4 7 10])
plot(X1(:,1),X1(:,2),'ro','LineWidth',3),grid
xlabel('x_i_1_j,  (i=1: j=1,...,n_1)  (i=2: j=1,...,n_2)','Fontsize',12); ylabel('x_i_2_j','Fontsize',12)
hold on
plot(X2(:,1),X2(:,2),'b*','LineWidth',3)
title({'Thus having "verified" the assumptions';'  ';'X_1 ~ N_2({\mu}_1,{\Sigma}_1)    (RED)';'X_2 ~ N_2({\mu}_2,{\Sigma}_2)   (BLUE)'},'Fontsize',14)
contour(y1,y2,F1,C1,'Color','r','LineWidth',1)
contour(z1,z2,F2,C2,'Color','b','LineWidth',1)
hold off
subplot(5,3,8:9)
axis off
line([0 1],[0.5 0.5],'Color','k','LineWidth',3)
subplot(5,3,[2 3 5 6])
axis off
Sp = ((n1-1)*S1 + (n2-1)*S2)/(n1+n2-2);
T = (n1+n2-2)*log(det(Sp)) - (n1-1)*log(det(S1)) - (n2-1)*log(det(S2));
correction_factor = 1 - ((2*p^2+3*p-1)/(6*(p+1)))*(1/(n1-1)+1/(n2-1)-1/(n1+n2-2));
test_statistic = correction_factor*T;
df = p*(p+1)/2;
p_value = 1-chi2cdf(test_statistic,df);
text(0,1,'Bartlett test of equal covariance matrices,  H_0: {\Sigma}_1 = {\Sigma}_2','Fontsize',20)
text(0,0.8,['S_1 = ' mat2str(S1,2), ',    n_1 = ' num2str(n1)],'Color','r','Fontsize',14)
text(0,0.7,['S_2 = ' mat2str(S2,2), ',    n_2 = ' num2str(n2)],'Color','b','Fontsize',14)
text(0,0.6,['S_p = [(n_1-1)S_1 + (n_2-1)S_2]/(n_1+n_2-2) = ' mat2str(Sp,2)],'Color','k','Fontsize',14)
text(0,0.4,'Test statistic, T = C*[(n_1+n_2-2)log|S_p| - (n_1-1)log|S_1| - (n_2-1)log|S_2|]  ~  {\chi}^2_p_(_p_+_1_)_/_2','Color','k','Fontsize',14)
text(0,0.28,'Finite population correction factor, C = 1 - ((2*p^2+3*p-1)/(6*(p+1)))*(1/(n_1-1)+1/(n_2-1)-1/(n_1+n_2-2))','Color','k','Fontsize',10)
text(0,0.1,['t = ' num2str(test_statistic,3)],'Color','k','Fontsize',14)
text(0,0,['p-value = P(T > t) = ' num2str(p_value,3) '   (where T ~ {\chi}^2_p_(_p_+_1_)_/_2)'],'Color','k','Fontsize',14)
text(0,-0.2,'CONCLUSION:  H_0 is (in spite of low p-value) ACCEPTED  =>  LINEAR Discriminant Analysis case','Fontsize',14)
%------------------------------------------------------------------------
% Hotelling T^2 test for significant separation of populations
%------------------------------------------------------------------------
subplot(5,3,[11 12 14 15])
axis off
test_statistic = ((n1+n2-p-1)/(p*(n1+n2-2)))*((x1_avg-x2_avg)*inv(Sp*(1/n1+1/n2))*(x1_avg-x2_avg)');
df1 = p;
df2 = n1+n2-p-1;
p_value = 1-fcdf(test_statistic,df1,df2);
text(0,1.1,'Hotelling T^2 test of significant separation of populations,  H_0: {\mu}_1 = {\mu}_2','Fontsize',20)
text(0,0.9,['x_1_,_a_v_g = ' mat2str(x1_avg,2), ',    n_1 = ' num2str(n1)],'Color','r','Fontsize',14)
text(0,0.8,['x_2_,_a_v_g = ' mat2str(x2_avg,2), ',    n_2 = ' num2str(n2)],'Color','b','Fontsize',14)
text(0,0.67,['S_p = ' mat2str(Sp,2)],'Color','k','Fontsize',14)
text(0,0.5,'Test statistic, T = [(n_1+n_2-p-1)/(p(n_1+n_2-2))] * [(x_1_,_a_v_g-x_2_,_a_v_g)^T[(1/n_1+1/n_2)S_p]^-^1(x_1_,_a_v_g-x_2_,_a_v_g)]  ~  F_p_,_n_1_+_n_2_-_p_-_1','Color','k','Fontsize',14)
text(0,0.33,['t = ' num2str(test_statistic,3)],'Color','k','Fontsize',14)
text(0,0.2,['p-value = P(T > t) = ' num2str(p_value,3) '   (where T ~ F_p_,_n_1_+_n_2_-_p_-_1)'],'Color','k','Fontsize',14)
text(0,0,'CONCLUSION:  H_0 is REJECTED  =>  Discriminant Analysis is relevant (significant separation of populations)','Fontsize',14)
%------------------------------------------------------------------------
% Linear discriminant function for equal priors and costs
%------------------------------------------------------------------------
figure('Name','Linear discriminant function for equal priors and costs')
subplot(4,7,1:3)
axis off
text(0,1.2,'Classification for two MVN populations with','Fontsize',22)
text(0,1,'   - equal misclassification costs, c(2|1) = c(1|2)','Fontsize',18)
text(0,0.8,'   - equal priors, p_1 = p_2','Fontsize',18)
text(0,0.5,'LINEAR discriminant decision rule for new observation x_0:','Fontsize',16)
text(0,0.25,'(x_1_,_a_v_g- x_2_,_a_v_g)^T S_p^-^1 [x_0 - 0.5(x_1_,_a_v_g+ x_2_,_a_v_g)]','Fontsize',16)
text(0.78,0.3,'> 0  {\rightarrow}  d_1','Fontsize',18)
text(0.78,0.15,'< 0  {\rightarrow}  d_2','Fontsize',18)
subplot(4,7,4:7:25)
axis off
line([0.5 0.5],[0 1],'Color','k','LineWidth',3)
annotation('arrow',[0.483125 0.559375],...
    [0.572979591836735 0.573979591836735],'HeadLength',20,'HeadWidth',20,...
    'LineWidth',6);
subplot(4,7,[8:10 15:17 22:24])
plot(X1(:,1),X1(:,2),'ro','LineWidth',3),grid,xlim([100 240]),ylim([300 750])
xlabel('x_i_1_j,  (i=1: j=1,...,n_1)  (i=2: j=1,...,n_2)','Fontsize',12); ylabel('x_i_2_j','Fontsize',12)
hold on
plot(X2(:,1),X2(:,2),'b*','LineWidth',3)
title('X_1~N_2({\mu}_1,{\Sigma}) (RED),   X_2~N_2({\mu}_2,{\Sigma}) (BLUE),   LDA-fct = 0 (BLACK)','Fontsize',12)
contour(y1,y2,F1,C1,'Color','r','LineWidth',1)
contour(z1,z2,F2,C2,'Color','b','LineWidth',1)
v1_min = min(min(y1),min(z1));
v1_max = max(max(y1),max(z1));
v1 = v1_min:(v1_max-v1_min)/100:v1_max;
v2_min = min(min(y2),min(z2));
v2_max = max(max(y2),max(z2));
v2 = v2_min:(v2_max-v2_min)/100:v2_max;
[V1,V2] = meshgrid(v1,v2);
F = (x1_avg-x2_avg)*inv(Sp)*[V1(:)' ; V2(:)'] - 0.5*(x1_avg-x2_avg)*inv(Sp)*(x1_avg+x2_avg)';
F = reshape(F,length(v2),length(v1));
contour(v1,v2,F,[eps eps],'Color','k','LineWidth',3)
coeffs_d1_LDA = [x1_avg*inv(Sp) -0.5*x1_avg*inv(Sp)*x1_avg'];
coeffs_d2_LDA = [x2_avg*inv(Sp) -0.5*x2_avg*inv(Sp)*x2_avg']; 
text(181,425,['d_1^L^D^A(x_1,x_2) = ' num2str(coeffs_d1_LDA(1),3) '*x_1 +'],'Fontsize',14)
text(212,398,[num2str(coeffs_d1_LDA(2),3) '*x_2 +'],'Fontsize',14)
text(212,380,num2str(coeffs_d1_LDA(3),3),'Fontsize',14)
text(110,675,['d_2^L^D^A(x_1,x_2) = ' num2str(coeffs_d2_LDA(1),3) '*x_1 +'],'Fontsize',14)
text(141,648,[num2str(coeffs_d2_LDA(2),3) '*x_2 +'],'Fontsize',14)
text(141,630,num2str(coeffs_d2_LDA(3),3),'Fontsize',14)
hold off
%------------------------------------------------------------------------
% Calculate confusion matrix and APER for training data
%------------------------------------------------------------------------
LDA_scores_pop1 = (x1_avg-x2_avg)*inv(Sp)*X1' - 0.5*(x1_avg-x2_avg)*inv(Sp)*(x1_avg+x2_avg)';
OK_classify_pop1 = find(LDA_scores_pop1 > 0);
n1C_APER = length(OK_classify_pop1);
not_OK_classify_pop1 = find(LDA_scores_pop1 < 0);
n1M_APER = length(not_OK_classify_pop1);
%----------------------------------
LDA_scores_pop2 = (x1_avg-x2_avg)*inv(Sp)*X2' - 0.5*(x1_avg-x2_avg)*inv(Sp)*(x1_avg+x2_avg)';
OK_classify_pop2 = find(LDA_scores_pop2 < 0);
n2C_APER = length(OK_classify_pop2);
not_OK_classify_pop2 = find(LDA_scores_pop2 > 0);
n2M_APER = length(not_OK_classify_pop2);
%----------------------------------
APER = (n1M_APER+n2M_APER)/(n1+n2);
%------------------------------------------------------------------------
% Calculate confusion matrix and AER_hat for training data using
% cross-validation (hold-out) procedure
%------------------------------------------------------------------------
n1C_AER = 0;
n1M_AER = 0;
for j1 = 1:n1
    X1_red = X1([1:(j1-1) j1+1:n1],:);
    x1_avg_red = mean(X1_red);
    S1_red = cov(X1_red);
    Sp_red = ((n1-2)*S1_red + (n2-1)*S2)/(n1+n2-3);
    LDA_score_heldout = (x1_avg_red-x2_avg)*inv(Sp_red)*X1(j1,:)' - 0.5*(x1_avg_red-x2_avg)*inv(Sp_red)*(x1_avg_red+x2_avg)';
    if LDA_score_heldout > 0
        n1C_AER = n1C_AER + 1;
    else
        n1M_AER = n1M_AER + 1;
    end
end
%----------------------------------
n2C_AER = 0;
n2M_AER = 0;
for j2 = 1:n2
    X2_red = X2([1:(j2-1) j2+1:n2],:);
    x2_avg_red = mean(X2_red);
    S2_red = cov(X2_red);
    Sp_red = ((n1-1)*S1 + (n2-2)*S2_red)/(n1+n2-3);
    LDA_score_heldout = (x1_avg-x2_avg_red)*inv(Sp_red)*X2(j2,:)' - 0.5*(x1_avg-x2_avg_red)*inv(Sp_red)*(x1_avg+x2_avg_red)';
    if LDA_score_heldout < 0
        n2C_AER = n2C_AER + 1;
    else
        n2M_AER = n2M_AER + 1;
    end
end
%----------------------------------
AER_hat = (n1M_AER+n2M_AER)/(n1+n2);
%------------------------------------------------------------------------
% Printout confusion tables and error rates
%------------------------------------------------------------------------
subplot(4,7,[5:7 12:14 19:21 26:28])
axis off
text(0,0.95,'Performance of classification','Fontsize',24)
text(0,0.9,'(empirical error rates)','Fontsize',24)
text(0,0.83,'------------------------------------------','Fontsize',22)
text(0,0.77,'Confusion matrix for training data:','Fontsize',20)
text(0.35,0.7,'classified as','Fontsize',18)
text(0.4,0.65,'d_1   d_2','Fontsize',18)
text(0.1,0.59,'true','Fontsize',18)
text(0.05,0.56,'population','Fontsize',18)
text(0.32,0.6,'{\pi}_1','Fontsize',18)
text(0.4,0.6,num2str(n1C_APER),'Fontsize',18)
text(0.5,0.6,num2str(n1M_APER),'Fontsize',18)
text(0.6,0.6,['(sum = n_1 = ' num2str(n1) ')'],'Fontsize',16)
text(0.32,0.55,'{\pi}_2','Fontsize',18)
text(0.42,0.55,num2str(n2M_APER),'Fontsize',18)
text(0.48,0.55,num2str(n2C_APER),'Fontsize',18)
text(0.6,0.55,['(sum = n_2 = ' num2str(n2) ')'],'Fontsize',16)
text(0,0.48,['Apparent Error Rate, APER = ' num2str(APER,2)],'Fontsize',18)
text(0,0.42,'------------------------------------------','Fontsize',22)
text(0,0.35,'Confusion matrix for training data','Fontsize',20)
text(0,0.31,'using cross-validation (hold-out):','Fontsize',20)
text(0.35,0.24,'classified as','Fontsize',18)
text(0.4,0.19,'d_1   d_2','Fontsize',18)
text(0.1,0.13,'true','Fontsize',18)
text(0.05,0.1,'population','Fontsize',18)
text(0.32,0.14,'{\pi}_1','Fontsize',18)
text(0.4,0.14,num2str(n1C_AER),'Fontsize',18)
text(0.5,0.14,num2str(n1M_AER),'Fontsize',18)
text(0.6,0.14,['(sum = n_1 = ' num2str(n1) ')'],'Fontsize',16)
text(0.32,0.09,'{\pi}_2','Fontsize',18)
text(0.42,0.09,num2str(n2M_AER),'Fontsize',18)
text(0.48,0.09,num2str(n2C_AER),'Fontsize',18)
text(0.6,0.09,['(sum = n_2 = ' num2str(n2) ')'],'Fontsize',16)
text(0,0.02,['Estimated Actual Error Rate, AER_h_a_t = ' num2str(AER_hat,2)],'Fontsize',18)
text(0,-0.04,'------------------------------------------','Fontsize',22)
%------------------------------------------------------------------------
% Classification of new observations with LDA
%------------------------------------------------------------------------
disp('---------------------------------------------------------------------')
disp('Classification of new observation 1 with LDA')
x0_1 = [170; 550]
d1_LDA_x0_1 = coeffs_d1_LDA*[x0_1; 1]
d2_LDA_x0_1 = coeffs_d2_LDA*[x0_1; 1]
disp('Conclusion: d2(x0_1) > d1(x0_1) => classify as member of population 2')
posterior_prob_x0_1_member_of_pop_1 = mvnpdf(x0_1,x1_avg',Sp)/(mvnpdf(x0_1,x1_avg',Sp) + mvnpdf(x0_1,x2_avg',Sp))
posterior_prob_x0_1_member_of_pop_2 = mvnpdf(x0_1,x2_avg',Sp)/(mvnpdf(x0_1,x1_avg',Sp) + mvnpdf(x0_1,x2_avg',Sp))
disp('Conclusion: P(pop 2 | x0_1) > P(pop 1 | x0_1) => classify as member of population 2')
disp('---------------------------------------------------------------------')
disp('Classification of new observation 2 with LDA')
x0_2 = [200; 600]
d1_LDA_x0_2 = coeffs_d1_LDA*[x0_2; 1]
d2_LDA_x0_2 = coeffs_d2_LDA*[x0_2; 1]
disp('Conclusion: d1(x0_2) > d2(x0_2) => classify as member of population 1')
posterior_prob_x0_2_member_of_pop_1 = mvnpdf(x0_2,x1_avg',Sp)/(mvnpdf(x0_2,x1_avg',Sp) + mvnpdf(x0_2,x2_avg',Sp))
posterior_prob_x0_2_member_of_pop_2 = mvnpdf(x0_2,x2_avg',Sp)/(mvnpdf(x0_2,x1_avg',Sp) + mvnpdf(x0_2,x2_avg',Sp))
disp('Conclusion: P(pop 1 | x0_2) > P(pop 2 | x0_2) => classify as member of population 1')
disp('---------------------------------------------------------------------')
disp('Classification of new observation 3 with LDA')
x0_3 = [120; 380]
d1_LDA_x0_3 = coeffs_d1_LDA*[x0_3; 1]
d2_LDA_x0_3 = coeffs_d2_LDA*[x0_3; 1]
disp('Conclusion: d1(x0_3) > d2(x0_3) => classify as member of population 1')
posterior_prob_x0_3_member_of_pop_1 = mvnpdf(x0_3,x1_avg',Sp)/(mvnpdf(x0_3,x1_avg',Sp) + mvnpdf(x0_3,x2_avg',Sp))
posterior_prob_x0_3_member_of_pop_2 = mvnpdf(x0_3,x2_avg',Sp)/(mvnpdf(x0_3,x1_avg',Sp) + mvnpdf(x0_3,x2_avg',Sp))
disp('Conclusion: P(pop 1 | x0_3) > P(pop 2 | x0_3) => classify as member of population 1')
disp('---------------------------------------------------------------------')
%-------------------------------------------------------------------------------

%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

%------------------------------------------------------------------------
% QDA: two populations, different mu, different SIGMA
%------------------------------------------------------------------------
% Bartlett test for equal covariance matrices for both populations
%------------------------------------------------------------------------
figure('Name','Bartlett and Hotelling tests --> QDA')
subplot(5,3,[4 7 10])
plot(X1(:,1),X1(:,2),'ro','LineWidth',3),grid
xlabel('x_i_1_j,  (i=1: j=1,...,n_1)  (i=2: j=1,...,n_2)','Fontsize',12); ylabel('x_i_2_j','Fontsize',12)
hold on
plot(X2(:,1),X2(:,2),'b*','LineWidth',3)
title({'Thus having "verified" the assumptions';'  ';'X_1 ~ N_2({\mu}_1,{\Sigma}_1)    (RED)';'X_2 ~ N_2({\mu}_2,{\Sigma}_2)   (BLUE)'},'Fontsize',14)
contour(y1,y2,F1,C1,'Color','r','LineWidth',1)
contour(z1,z2,F2,C2,'Color','b','LineWidth',1)
hold off
subplot(5,3,8:9)
axis off
line([0 1],[0.5 0.5],'Color','k','LineWidth',3)
subplot(5,3,[2 3 5 6])
axis off
Sp = ((n1-1)*S1 + (n2-1)*S2)/(n1+n2-2);
T = (n1+n2-2)*log(det(Sp)) - (n1-1)*log(det(S1)) - (n2-1)*log(det(S2));
correction_factor = 1 - ((2*p^2+3*p-1)/(6*(p+1)))*(1/(n1-1)+1/(n2-1)-1/(n1+n2-2));
test_statistic = correction_factor*T;
df = p*(p+1)/2;
p_value = 1-chi2cdf(test_statistic,df);
text(0,1,'Bartlett test of equal covariance matrices,  H_0: {\Sigma}_1 = {\Sigma}_2','Fontsize',20)
text(0,0.8,['S_1 = ' mat2str(S1,2), ',    n_1 = ' num2str(n1)],'Color','r','Fontsize',14)
text(0,0.7,['S_2 = ' mat2str(S2,2), ',    n_2 = ' num2str(n2)],'Color','b','Fontsize',14)
text(0,0.6,['S_p = [(n_1-1)S_1 + (n_2-1)S_2]/(n_1+n_2-2) = ' mat2str(Sp,2)],'Color','k','Fontsize',14)
text(0,0.4,'Test statistic, T = C*[(n_1+n_2-2)log|S_p| - (n_1-1)log|S_1| - (n_2-1)log|S_2|]  ~  {\chi}^2_p_(_p_+_1_)_/_2','Color','k','Fontsize',14)
text(0,0.28,'Finite population correction factor, C = 1 - ((2*p^2+3*p-1)/(6*(p+1)))*(1/(n_1-1)+1/(n_2-1)-1/(n_1+n_2-2))','Color','k','Fontsize',10)
text(0,0.1,['t = ' num2str(test_statistic,3)],'Color','k','Fontsize',14)
text(0,0,['p-value = P(T > t) = ' num2str(p_value,3) '   (where T ~ {\chi}^2_p_(_p_+_1_)_/_2)'],'Color','k','Fontsize',14)
text(0,-0.2,'CONCLUSION:  H_0 is REJECTED  =>  QUADRATIC Discriminant Analysis case','Fontsize',14)
%------------------------------------------------------------------------
% Hotelling T^2 test for significant separation of populations
%------------------------------------------------------------------------
subplot(5,3,[11 12 14 15])
axis off
test_statistic = (x1_avg-x2_avg)*inv(S1/n1+S2/n2)*(x1_avg-x2_avg)';
df = p;
p_value = 1-chi2cdf(test_statistic,df);
text(0,1.1,'Hotelling T^2 test of significant separation of populations,  H_0: {\mu}_1 = {\mu}_2','Fontsize',20)
text(0,0.9,['x_1_,_a_v_g = ' mat2str(x1_avg,2), ',    n_1 = ' num2str(n1)],'Color','r','Fontsize',14)
text(0,0.8,['x_2_,_a_v_g = ' mat2str(x2_avg,2), ',    n_2 = ' num2str(n2)],'Color','b','Fontsize',14)
text(0,0.67,['S_1 = ' mat2str(S1,2)],'Color','r','Fontsize',14)
text(0,0.57,['S_2 = ' mat2str(S2,2)],'Color','b','Fontsize',14)
text(0,0.4,'Test statistic, T = [(x_1_,_a_v_g-x_2_,_a_v_g)^T(S_1/n_1+S_2/n_2)^-^1(x_1_,_a_v_g-x_2_,_a_v_g)]  ~  {\chi}^2_p','Color','k','Fontsize',14)
text(0,0.23,['t = ' num2str(test_statistic,3)],'Color','k','Fontsize',14)
text(0,0.1,['p-value = P(T > t) = ' num2str(p_value,3) '   (where T ~ {\chi}^2_p)'],'Color','k','Fontsize',14)
text(0,-0.1,'CONCLUSION:  H_0 is REJECTED  =>  Discriminant Analysis is relevant (significant separation of populations)','Fontsize',14)
%------------------------------------------------------------------------
% Quadratic discriminant function for equal priors and costs
%------------------------------------------------------------------------
figure('Name','Quadratic discriminant function for equal priors and costs')
subplot(4,7,1:3)
axis off
text(-0.25,1.2,'Classification for two MVN populations with equal misclassification costs, c(2|1) = c(1|2) and equal priors, p_1 = p_2','Fontsize',22)
text(-0.1,0.8,'QUADRATIC discriminant decision rule for new observation x_0:','Fontsize',16)
text(-0.1,0.5,'-0.5*x_0^T(S_1^-^1- S_2^-^1)x_0 + (x_1_,_a_v_g^TS_1^-^1- x_2_,_a_v_g^TS_2^-^1)x_0 - k','Fontsize',16)
text(-0.1,0.18,' where  k = 0.5*log(|S_1|/|S_2|) + 0.5*(x_1_,_a_v_g^TS_1^-^1x_1_,_a_v_g- x_2_,_a_v_g^TS_2^-^1x_2_,_a_v_g)','Fontsize',14)
text(0.76,0.55,'> 0  {\rightarrow}  d_1','Fontsize',18)
text(0.76,0.4,'< 0  {\rightarrow}  d_2','Fontsize',18)
subplot(4,7,4:7:25)
axis off
line([0.5 0.5],[0 1],'Color','k','LineWidth',3)
annotation('arrow',[0.483125 0.559375],...
    [0.572979591836735 0.573979591836735],'HeadLength',20,'HeadWidth',20,...
    'LineWidth',6);
subplot(4,7,[8:10 15:17 22:24])
plot(X1(:,1),X1(:,2),'ro','LineWidth',3),grid,xlim([100 240]),ylim([300 750])
xlabel('x_i_1_j,  (i=1: j=1,...,n_1)  (i=2: j=1,...,n_2)','Fontsize',12); ylabel('x_i_2_j','Fontsize',12)
hold on
plot(X2(:,1),X2(:,2),'b*','LineWidth',3)
title('X_1~N_2({\mu}_1,{\Sigma}) (RED),   X_2~N_2({\mu}_2,{\Sigma}) (BLUE),   QDA-fct = 0 (BLACK)','Fontsize',12)
contour(y1,y2,F1,C1,'Color','r','LineWidth',1)
contour(z1,z2,F2,C2,'Color','b','LineWidth',1)
v1_min = min(min(y1),min(z1));
v1_max = max(max(y1),max(z1));
v1 = v1_min:(v1_max-v1_min)/100:v1_max;
v2_min = min(min(y2),min(z2));
v2_max = max(max(y2),max(z2));
v2 = v2_min:(v2_max-v2_min)/100:v2_max;
[V1,V2] = meshgrid(v1,v2);
K = 0.5*log(det(S1)/det(S2)) + 0.5*(x1_avg*inv(S1)*x1_avg' - x2_avg*inv(S2)*x2_avg');
V1_vector = V1(:);
V2_vector = V2(:);
F_QDA = zeros(1,length(V1_vector));
for q = 1:length(V1_vector)
    F_QDA(q) = -0.5*[V1(q) V2(q)]*(inv(S1) - inv(S2))*[V1(q); V2(q)] + (x1_avg*inv(S1) - x2_avg*inv(S2))*[V1(q); V2(q)] - K;
end
F_QDA = reshape(F_QDA,length(v2),length(v1));
contour(v1,v2,F_QDA,[eps eps],'Color','k','LineWidth',3)
quad_coeffs_d1_QDA = -0.5*inv(S1);
lin_coeffs_d1_QDA = x1_avg*inv(S1);
const_coeff_d1_QDA = -0.5*log(det(S1)) - 0.5*x1_avg*inv(S1)*x1_avg';
quad_coeffs_d2_QDA = -0.5*inv(S2);
lin_coeffs_d2_QDA = x2_avg*inv(S2);
const_coeff_d2_QDA = -0.5*log(det(S2)) - 0.5*x2_avg*inv(S2)*x2_avg';
text(170,425,['d_1^Q^D^A(x_1,x_2) = ' num2str(quad_coeffs_d1_QDA(1,1),3) '*x_1^2 +'],'Fontsize',14)
text(202,400,[num2str(quad_coeffs_d1_QDA(2,2),3) '*x_2^2 +'],'Fontsize',14)
text(202,375,[num2str(2*quad_coeffs_d1_QDA(1,2),3) '*x_1*x_2 +'],'Fontsize',14)
text(202,352,[num2str(lin_coeffs_d1_QDA(1),3) '*x_1 +'],'Fontsize',14)
text(202,329,[num2str(lin_coeffs_d1_QDA(2),3) '*x_2 +'],'Fontsize',14)
text(202,309,num2str(const_coeff_d1_QDA,3),'Fontsize',14)
text(105,725,['d_2^Q^D^A(x_1,x_2) = ' num2str(quad_coeffs_d2_QDA(1,1),3) '*x_1^2 +'],'Fontsize',14)
text(137,700,[num2str(quad_coeffs_d2_QDA(2,2),3) '*x_2^2 +'],'Fontsize',14)
text(137,675,[num2str(2*quad_coeffs_d2_QDA(1,2),3) '*x_1*x_2 +'],'Fontsize',14)
text(137,652,[num2str(lin_coeffs_d2_QDA(1),3) '*x_1 +'],'Fontsize',14)
text(137,629,[num2str(lin_coeffs_d2_QDA(2),3) '*x_2 +'],'Fontsize',14)
text(137,609,num2str(const_coeff_d2_QDA,3),'Fontsize',14)
hold off
%------------------------------------------------------------------------
% Calculate confusion matrix and APER for training data
%------------------------------------------------------------------------
QDA_scores_pop1 = zeros(1,n1);
for q1 = 1:n1
    QDA_scores_pop1(q1) = -0.5*X1(q1,:)*(inv(S1) - inv(S2))*X1(q1,:)' + (x1_avg*inv(S1) - x2_avg*inv(S2))*X1(q1,:)' - K;
end    
OK_classify_pop1 = find(QDA_scores_pop1 > 0);
n1C_APER = length(OK_classify_pop1);
not_OK_classify_pop1 = find(QDA_scores_pop1 < 0);
n1M_APER = length(not_OK_classify_pop1);
%----------------------------------
QDA_scores_pop2 = zeros(1,n2);
for q2 = 1:n2
    QDA_scores_pop2(q2) = -0.5*X2(q2,:)*(inv(S1) - inv(S2))*X2(q2,:)' + (x1_avg*inv(S1) - x2_avg*inv(S2))*X2(q2,:)' - K;
end  
OK_classify_pop2 = find(QDA_scores_pop2 < 0);
n2C_APER = length(OK_classify_pop2);
not_OK_classify_pop2 = find(QDA_scores_pop2 > 0);
n2M_APER = length(not_OK_classify_pop2);
%----------------------------------
APER = (n1M_APER+n2M_APER)/(n1+n2);
%------------------------------------------------------------------------
% Calculate confusion matrix and AER_hat for training data using
% cross-validation (hold-out) procedure
%------------------------------------------------------------------------
n1C_AER = 0;
n1M_AER = 0;
for q1 = 1:n1
    X1_red = X1([1:(q1-1) q1+1:n1],:);
    x1_avg_red = mean(X1_red);
    S1_red = cov(X1_red);
    K_red = 0.5*log(det(S1_red)/det(S2)) + 0.5*(x1_avg_red*inv(S1_red)*x1_avg_red' - x2_avg*inv(S2)*x2_avg');
    QDA_score_heldout = -0.5*X1(q1,:)*(inv(S1_red) - inv(S2))*X1(q1,:)' + (x1_avg_red*inv(S1_red) - x2_avg*inv(S2))*X1(q1,:)' - K_red;
    if QDA_score_heldout > 0
        n1C_AER = n1C_AER + 1;
    else
        n1M_AER = n1M_AER + 1;
    end
end
%----------------------------------
n2C_AER = 0;
n2M_AER = 0;
for q2 = 1:n2
    X2_red = X2([1:(q2-1) q2+1:n2],:);
    x2_avg_red = mean(X2_red);
    S2_red = cov(X2_red);
    K_red = 0.5*log(det(S1)/det(S2_red)) + 0.5*(x1_avg*inv(S1)*x1_avg' - x2_avg_red*inv(S2_red)*x2_avg_red');
    QDA_score_heldout = -0.5*X2(q2,:)*(inv(S1) - inv(S2_red))*X2(q2,:)' + (x1_avg*inv(S1) - x2_avg_red*inv(S2_red))*X2(q2,:)' - K_red;
    if QDA_score_heldout < 0
        n2C_AER = n2C_AER + 1;
    else
        n2M_AER = n2M_AER + 1;
    end
end
%----------------------------------
AER_hat = (n1M_AER+n2M_AER)/(n1+n2);
%------------------------------------------------------------------------
% Printout confusion tables and error rates
%------------------------------------------------------------------------
subplot(4,7,[5:7 12:14 19:21 26:28])
axis off
text(0,0.95,'Performance of classification','Fontsize',24)
text(0,0.9,'(empirical error rates)','Fontsize',24)
text(0,0.83,'------------------------------------------','Fontsize',22)
text(0,0.77,'Confusion matrix for training data:','Fontsize',20)
text(0.35,0.7,'classified as','Fontsize',18)
text(0.4,0.65,'d_1   d_2','Fontsize',18)
text(0.1,0.59,'true','Fontsize',18)
text(0.05,0.56,'population','Fontsize',18)
text(0.32,0.6,'{\pi}_1','Fontsize',18)
text(0.4,0.6,num2str(n1C_APER),'Fontsize',18)
text(0.5,0.6,num2str(n1M_APER),'Fontsize',18)
text(0.6,0.6,['(sum = n_1 = ' num2str(n1) ')'],'Fontsize',16)
text(0.32,0.55,'{\pi}_2','Fontsize',18)
text(0.42,0.55,num2str(n2M_APER),'Fontsize',18)
text(0.48,0.55,num2str(n2C_APER),'Fontsize',18)
text(0.6,0.55,['(sum = n_2 = ' num2str(n2) ')'],'Fontsize',16)
text(0,0.48,['Apparent Error Rate, APER = ' num2str(APER,2)],'Fontsize',18)
text(0,0.42,'------------------------------------------','Fontsize',22)
text(0,0.35,'Confusion matrix for training data','Fontsize',20)
text(0,0.31,'using cross-validation (hold-out):','Fontsize',20)
text(0.35,0.24,'classified as','Fontsize',18)
text(0.4,0.19,'d_1   d_2','Fontsize',18)
text(0.1,0.13,'true','Fontsize',18)
text(0.05,0.1,'population','Fontsize',18)
text(0.32,0.14,'{\pi}_1','Fontsize',18)
text(0.4,0.14,num2str(n1C_AER),'Fontsize',18)
text(0.5,0.14,num2str(n1M_AER),'Fontsize',18)
text(0.6,0.14,['(sum = n_1 = ' num2str(n1) ')'],'Fontsize',16)
text(0.32,0.09,'{\pi}_2','Fontsize',18)
text(0.42,0.09,num2str(n2M_AER),'Fontsize',18)
text(0.48,0.09,num2str(n2C_AER),'Fontsize',18)
text(0.6,0.09,['(sum = n_2 = ' num2str(n2) ')'],'Fontsize',16)
text(0,0.02,['Estimated Actual Error Rate, AER_h_a_t = ' num2str(AER_hat,2)],'Fontsize',18)
text(0,-0.04,'------------------------------------------','Fontsize',22)
%------------------------------------------------------------------------
% Classification of new observations with QDA
%------------------------------------------------------------------------
disp('    ')
disp('---------------------------------------------------------------------')
disp('Classification of new observation 1 with QDA')
x0_1
d1_QDA_x0_1 = x0_1'*quad_coeffs_d1_QDA*x0_1 + lin_coeffs_d1_QDA*x0_1 + const_coeff_d1_QDA
d2_QDA_x0_1 = x0_1'*quad_coeffs_d2_QDA*x0_1 + lin_coeffs_d2_QDA*x0_1 + const_coeff_d2_QDA
disp('Conclusion: d2(x0_1) > d1(x0_1) => classify as member of population 2')
posterior_prob_x0_1_member_of_pop_1 = mvnpdf(x0_1,x1_avg',S1)/(mvnpdf(x0_1,x1_avg',S1) + mvnpdf(x0_1,x2_avg',S2))
posterior_prob_x0_1_member_of_pop_2 = mvnpdf(x0_1,x2_avg',S2)/(mvnpdf(x0_1,x1_avg',S1) + mvnpdf(x0_1,x2_avg',S2))
disp('Conclusion: P(pop 2 | x0_1) > P(pop 1 | x0_1) => classify as member of population 2 (same as LDA)')
disp('---------------------------------------------------------------------')
disp('Classification of new observation 2 with QDA')
x0_2
d1_QDA_x0_2 = x0_2'*quad_coeffs_d1_QDA*x0_2 + lin_coeffs_d1_QDA*x0_2 + const_coeff_d1_QDA
d2_QDA_x0_2 = x0_2'*quad_coeffs_d2_QDA*x0_2 + lin_coeffs_d2_QDA*x0_2 + const_coeff_d2_QDA
disp('Conclusion: d1(x0_2) > d2(x0_2) => classify as member of population 1')
posterior_prob_x0_2_member_of_pop_1 = mvnpdf(x0_2,x1_avg',S1)/(mvnpdf(x0_2,x1_avg',S1) + mvnpdf(x0_2,x2_avg',S2))
posterior_prob_x0_2_member_of_pop_2 = mvnpdf(x0_2,x2_avg',S2)/(mvnpdf(x0_2,x1_avg',S1) + mvnpdf(x0_2,x2_avg',S2))
disp('Conclusion: P(pop 1 | x0_2) > P(pop 2 | x0_2) => classify as member of population 1 (same as LDA) ')
disp('---------------------------------------------------------------------')
disp('Classification of new observation 3 with QDA')
x0_3
d1_QDA_x0_3 = x0_3'*quad_coeffs_d1_QDA*x0_3 + lin_coeffs_d1_QDA*x0_3 + const_coeff_d1_QDA
d2_QDA_x0_3 = x0_3'*quad_coeffs_d2_QDA*x0_3 + lin_coeffs_d2_QDA*x0_3 + const_coeff_d2_QDA
disp('Conclusion: d2(x0_3) > d1(x0_3) => classify as member of population 2')
posterior_prob_x0_3_member_of_pop_1 = mvnpdf(x0_3,x1_avg',S1)/(mvnpdf(x0_3,x1_avg',S1) + mvnpdf(x0_3,x2_avg',S2))
posterior_prob_x0_3_member_of_pop_2 = mvnpdf(x0_3,x2_avg',S2)/(mvnpdf(x0_3,x1_avg',S1) + mvnpdf(x0_3,x2_avg',S2))
disp('Conclusion: P(pop 2 | x0_3) > P(pop 1 | x0_3) => classify as member of population 2 (different from LDA)')
disp('---------------------------------------------------------------------')
%-------------------------------------------------------------------------------------------------------