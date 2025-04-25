clc;clear;close all
format compact

disp('--------------------------------------')
X1 = [0.2; 0.8];
X2 = [0.25; 0.7];
X3 = [0.4; 0.7];
X4 = [0.6; 0.35];
X5 = [0.8; 0.4];
X_cluster_1 = [X1'; X2'; X3']
X_cluster_2 = [X4'; X5']
disp('--------------------------------------')
centroid_1 = mean(X_cluster_1)
centroid_2 = mean(X_cluster_2)
disp('--------------------------------------')
d14 = norm(X1-X4)
d15 = norm(X1-X5)
d24 = norm(X2-X4)
d25 = norm(X2-X5)
d34 = norm(X3-X4)
d35 = norm(X3-X5)
disp('--------------------------------------')
d12_single_linkage = min([d14 d15 d24 d25 d34 d35])
d12_complete_linkage = max([d14 d15 d24 d25 d34 d35])
d12_average_linkage = mean([d14 d15 d24 d25 d34 d35])
d12_centroid_linkage = norm(centroid_1 - centroid_2)
disp('--------------------------------------')

figure
subplot(2,4,1:4)
axis off
text(0,0.9,'Inter-cluster distance metrics: Linkage methods','Fontsize',28)
text(0,0.78,'(using Euclidean norm as point-to-point distance metric)','Fontsize',20)
text(0,0.5,'Single linkage','Fontsize',20)
text(0,0.4,'(nearest neighbor)','Fontsize',18)
text(0,0.17,'d_1_2 = min (d_i_j), i{\in}C_1, j{\in}C_2','Fontsize',16)
text(0,-0.05,'example:','Fontsize',16)
text(0,-0.15,['d_1_2 = ' num2str(d12_single_linkage,4)],'Fontsize',16)
text(0.27,0.5,'Complete linkage','Fontsize',20)
text(0.27,0.4,'(farthest neighbor)','Fontsize',18)
text(0.27,0.17,'d_1_2 = max (d_i_j), i{\in}C_1, j{\in}C_2','Fontsize',16)
text(0.27,-0.05,'example:','Fontsize',16)
text(0.27,-0.15,['d_1_2 = ' num2str(d12_complete_linkage,4)],'Fontsize',16)
text(0.54,0.5,'Average linkage','Fontsize',20)
text(0.54,0.4,'(average distance)','Fontsize',18)
text(0.54,0.17,'d_1_2 = ({\Sigma}_i{\Sigma}_j d_i_j)/N_1N_2, i{\in}C_1, j{\in}C_2','Fontsize',16)
text(0.54,-0.05,'example:','Fontsize',16)
text(0.54,-0.15,['d_1_2 = ' num2str(d12_average_linkage,4)],'Fontsize',16)
text(0.81,0.5,'Centroid linkage','Fontsize',20)
text(0.81,0.4,'(centroid distance)','Fontsize',18)
text(0.81,0.17,'d_1_2 = d (X_i,X_j), i{\in}C_1, j{\in}C_2','Fontsize',16)
text(0.813,0.27,'             _  _','Fontsize',16),warning off
text(0.81,-0.05,'example:','Fontsize',16)
text(0.81,-0.15,['d_1_2 = ' num2str(d12_centroid_linkage,4)],'Fontsize',16)

subplot(2,4,5)
plot(X_cluster_1(:,1),X_cluster_1(:,2),'bo','LineWidth',5),xlim([0 1]),ylim([0 1]),grid,axis square
hold on
plot(X_cluster_2(:,1),X_cluster_2(:,2),'ro','LineWidth',5)
text(0.22,0.84,'X1')
text(0.16,0.7,'X2')
text(0.38,0.76,'X3')
text(0.57,0.29,'X4')
text(0.82,0.36,'X5')
annotation('ellipse',[0.14 0.302295918367347 0.0737500000000002 0.119897959183673]);
annotation('ellipse',[0.210375 0.17984693877551 0.0615 0.105867346938776]);
text(0.5,0.9,'Cluster 1','Fontsize',12)
text(0.5,0.1,'Cluster 2','Fontsize',12)
annotation('doublearrow',[0.19375 0.220625],[0.339285714285714 0.241071428571429],'LineWidth',3);
hold off

subplot(2,4,6)
plot(X_cluster_1(:,1),X_cluster_1(:,2),'bo','LineWidth',5),xlim([0 1]),,ylim([0 1]),grid,axis square
hold on
plot(X_cluster_2(:,1),X_cluster_2(:,2),'ro','LineWidth',5)
text(0.22,0.84,'X1')
text(0.16,0.7,'X2')
text(0.38,0.76,'X3')
text(0.57,0.29,'X4')
text(0.82,0.36,'X5')
annotation('ellipse',[0.346875000000001 0.302704081632653 0.0737500000000002 0.119897959183673]);
annotation('ellipse',[0.414125 0.181530612244897 0.0615 0.105867346938776]);
text(0.5,0.9,'Cluster 1','Fontsize',12)
text(0.5,0.1,'Cluster 2','Fontsize',12)
annotation('doublearrow',[0.37 0.459375],[0.371448979591837 0.252551020408163],'LineWidth',3);
hold off

subplot(2,4,7)
plot(X_cluster_1(:,1),X_cluster_1(:,2),'bo','LineWidth',5),xlim([0 1]),,ylim([0 1]),grid,axis square
hold on
plot(X_cluster_2(:,1),X_cluster_2(:,2),'ro','LineWidth',5)
text(0.22,0.84,'X1')
text(0.16,0.7,'X2')
text(0.38,0.76,'X3')
text(0.57,0.29,'X4')
text(0.82,0.36,'X5')
annotation('ellipse',[0.555000000000002 0.297602040816325 0.0737500000000002 0.119897959183673]);
annotation('ellipse',[0.621625000000001 0.178979591836734 0.0615 0.105867346938776]);
text(0.5,0.9,'Cluster 1','Fontsize',12)
text(0.5,0.1,'Cluster 2','Fontsize',12)
annotation('doublearrow',[0.5825 0.634375],[0.337010204081633 0.232142857142857],'LineWidth',3);
annotation('doublearrow',[0.575 0.635625],[0.370173469387755 0.235969387755102],'LineWidth',3);
annotation('doublearrow',[0.606875 0.633125],[0.334459183673469 0.241071428571429],'LineWidth',3);
annotation('doublearrow',[0.584375 0.66375],[0.335734693877551 0.246173469387755],'LineWidth',3);
annotation('doublearrow',[0.576875 0.665625],[0.370173469387755 0.247448979591837],'LineWidth',3);
annotation('doublearrow',[0.6075 0.664375],[0.337010204081633 0.251275510204082],'LineWidth',3);
hold off

subplot(2,4,8)
plot(X_cluster_1(:,1),X_cluster_1(:,2),'bo','LineWidth',5),xlim([0 1]),,ylim([0 1]),grid,axis square
hold on
plot(X_cluster_2(:,1),X_cluster_2(:,2),'ro','LineWidth',5)
plot(centroid_1(1),centroid_1(2),'gx','LineWidth',10)
plot(centroid_2(1),centroid_2(2),'gx','LineWidth',10)
text(0.22,0.84,'X1')
text(0.16,0.7,'X2')
text(0.38,0.76,'X3')
text(0.57,0.29,'X4')
text(0.82,0.36,'X5')
annotation('ellipse',[0.763125 0.295051020408163 0.0737500000000002 0.119897959183673]);
annotation('ellipse',[0.831 0.177704081632652 0.0615 0.105867346938776]);
text(0.5,0.9,'Cluster 1','Fontsize',12)
text(0.5,0.1,'Cluster 2','Fontsize',12)
annotation('doublearrow',[0.794375 0.855],[0.349765306122449 0.247448979591837],'LineWidth',3);
hold off
%-------------------------------------------------------------------------------------------------------