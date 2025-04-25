clear;close all;clc;
format compact

load('M:\Statistik\MULTISTAT\Lektion 12  Cluster Analysis\Course material\dataset_mds_demo');
load('M:\Statistik\MULTISTAT\Lektion 12  Cluster Analysis\Course material\countriesID_mds_demo');

disp('---------------------------------------------------------------------')
disp('Descriptive statistics estimators for female speed-data')
disp('---------------------------------------------------------------------')
[n p] = size(X_time)
X_speed(:,1) = 100./X_time(:,1);
X_speed(:,2) = 200./X_time(:,2);
X_speed(:,3) = 400./X_time(:,3);
X_speed(:,4) = 800./(60*X_time(:,4));
X_speed(:,5) = 1500./(60*X_time(:,5));
X_speed(:,6) = 3000./(60*X_time(:,6));
X_speed(:,7) = 42195./(60*X_time(:,7));
Xmean_speed = mean(X_speed,1);
S_speed = cov(X_speed)
R_speed = corr(X_speed)

disp('---------------------------------------------------------------------')
disp('Multidimensional Scaling plot for female speed-data')
disp('---------------------------------------------------------------------')
D = pdist(X_speed,'euclidean');
[Y eigvals] = cmdscale(D);
biggest_eigvalues = [eigvals(1:10) eigvals(1:10)/max(abs(eigvals))]
figure('Name','Multidimensional Scaling plot for female speed-data')
subplot(5,1,1)
axis off
text(0.03,1,'Classical metric multidimensional scaling, example','Fontsize',28)
text(0.03,0.7,'Approximately distance-preserving plot of p=7 dimensional data in q=2 dimensions','Fontsize',24)
text(0.03,0.2,'Example data:','Fontsize',18)
text(0.03,0,'Female national track records in the 7 disciplines: 100m, 200m, 400m, 800m, 1500m, 3000m, Marathon','Fontsize',18)
text(0.03,-0.2,'(all variables in unit [m/s])','Fontsize',18)
subplot(5,1,2:5)
plot(Y(:,1),Y(:,2),'.'),for j = 1:54 text(Y(j,1)+0.02,Y(j,2),countriesIDmdsdemo(j,1),'Fontsize',8); end,xlabel('MDS dim #1 scores','Fontsize',14),ylabel('MDS dim #2 scores','Fontsize',14)
    
disp('-------------------------------------------------------')
disp('PCA on S for female speed-data')
disp('-------------------------------------------------------')
[E_S_speed LAMBDA_S_speed E_S_speed] = svd(S_speed);
E_S_speed = sign(E_S_speed(1,1))*E_S_speed;
LAMBDA_S_speed,E_S_speed
PCA_1_S_speed_loadings_X_variables = (E_S_speed(:,1))'
PCA_1_S_speed_scores_X_variables = (X_speed - ones(n,1)*mean(X_speed))*PCA_1_S_speed_loadings_X_variables';
PCA_2_S_speed_loadings_X_variables = (E_S_speed(:,2))'
PCA_2_S_speed_scores_X_variables = (X_speed - ones(n,1)*mean(X_speed))*PCA_2_S_speed_loadings_X_variables';

figure('Name','Scatterplots for female speed-data PCA 1 and 2 scores')
subplot(5,3,1:3)
axis off
text(0.05,1.2,'Equivalence between','Fontsize',24)
text(0.05,0.9,'Classical metric Multidimensional Scaling (MDS) and','Fontsize',24)
text(0.05,0.6,'Principal Component Analysis (PCA)','Fontsize',24)
text(0.03,-0.05,'MDS scores for q=2','Fontsize',20,'color','b')
text(0.4,-0.05,'PC_1 and PC_2 scores','Fontsize',20,'color','r')
text(0.75,-0.05,'- PC_1 and - PC_2 scores','Fontsize',20,'color','m')
text(0.78,-0.27,'(reflection of both axes)','Fontsize',14,'color','m')

subplot(5,3,5:3:14)
plot(PCA_1_S_speed_scores_X_variables,PCA_2_S_speed_scores_X_variables,'r.'),xlim([-3 1.5]),ylim([-1 0.5])
xlabel('1^s^t Principal Component scores','Fontsize',14,'color','r')
ylabel('2^n^d Principal Component scores','Fontsize',14,'color','r')
hold on
for j = 1:54 text(PCA_1_S_speed_scores_X_variables(j)+0.02,PCA_2_S_speed_scores_X_variables(j),countriesIDmdsdemo(j,1),'Fontsize',8); end
hold off
subplot(5,3,6:3:15)
plot(-PCA_1_S_speed_scores_X_variables,-PCA_2_S_speed_scores_X_variables,'m.'),,xlim([-1.5 3]),ylim([-0.5 1])
xlabel('- 1^s^t Principal Component scores','Fontsize',14,'color','m')
ylabel('- 2^n^d Principal Component scores','Fontsize',14,'color','m')
hold on
for j = 1:54 text(-PCA_1_S_speed_scores_X_variables(j)+0.02,-PCA_2_S_speed_scores_X_variables(j),countriesIDmdsdemo(j,1),'Fontsize',8); end
hold off
subplot(5,3,4:3:13)
plot(Y(:,1),Y(:,2),'.'),for j = 1:54 text(Y(j,1)+0.02,Y(j,2),countriesIDmdsdemo(j,1),'Fontsize',8); end,...
xlabel('MDS dim #1 scores','Fontsize',14,'color','b'),ylabel('MDS dim #2 scores','Fontsize',14,'color','b'),xlim([-1.5 3]),ylim([-0.5 1])
   
%----------------------------------------------------

figure('Name','Bi-plot for female speed-data PCA 1 and 2 scores and components') 
subplot(5,1,1)
axis off
text(0.05,0.75,'Biplot:','Fontsize',28)
text(0.05,0.4,'Joined plot of PC scores for observations and PC loadings on variables','Fontsize',20)
text(0.05,0.1,'(note that "biplot" scales observation scores and also determines sign of loadings for "nice plotting")','Fontsize',14)
subplot(5,1,2:5)
names = char('100','200','400','800','1500','3000','Marathon');
biplot(E_S_speed(:,1:2),'scores',[PCA_1_S_speed_scores_X_variables PCA_2_S_speed_scores_X_variables],'varlabels',names)
xlabel('PC_1 scores for observations and loadings on variables','Fontsize',14,'color','k')
ylabel('PC_2 scores for observations and loadings on variables','Fontsize',14,'color','k')
disp('-------------------------------------------------------')
%------------------------------------------------------------------------------------------------------------------------------

