clear;close all;clc;
format compact

disp('-------------------------------------------------------------------------------------------------------')
X = [1.1 9.1; 1.5 7.9; 2.1 7.5; 2.9 9.2; 0.5 0.9; 2.4 4.7; 3.1 5.4; 3.8 3.5; 7.5 3.2; 8.8 3.9]
[n p] = size(X)
disp('-------------------------------------------------------------------------------------------------------')
D = squareform(pdist(X,'euclidean'))
disp('-------------------------------------------------------------------------------------------------------')
Tree = linkage(D,'single');
disp('assigments with 3 clusters (single linkage):')
assignments_3_clusters = cluster(Tree,'maxclust',3);
cluster_l = find(assignments_3_clusters == 1)'
cluster_2 = find(assignments_3_clusters == 2)'
cluster_3 = find(assignments_3_clusters == 3)'
disp('-------------------------------------------------------------------------------------------------------')
disp('assigments with 4 clusters (single linkage):')
assignments_4_clusters = cluster(Tree,'maxclust',4);
cluster_l = find(assignments_4_clusters == 1)'
cluster_2 = find(assignments_4_clusters == 2)'
cluster_3 = find(assignments_4_clusters == 3)'
cluster_4 = find(assignments_4_clusters == 4)'
disp('-------------------------------------------------------------------------------------------------------')

figure
subplot(5,2,1:2)
axis off
text(0,1,'Hierarchical agglomerative clustering example, using','Fontsize',20)
text(0.53,1,'-  Euclidean point-to-point distance metric','Fontsize',18)
text(0.53,0.75,'-  Single linkage (nearest neighbor) inter-cluster distance metric','Fontsize',18)

subplot(5,2,3:2:9)
plot(X(:,1),X(:,2),'ro','LineWidth',5),xlim([0 10]),ylim([0 10]),grid,axis square,xlabel('x_1','Fontsize',14),ylabel('x_2','Fontsize',14)
title({'10 bivariate observations and'; '2 possible clusterings based on dendrogram cuts'},'Fontsize',14)
hold on
for j = 1:n
    text(X(j,1)+0.2,X(j,2),num2str(j),'Fontsize',14)
end
annotation('ellipse',[0.155375 0.537397959183669 0.104 0.211734693877547],'LineStyle','--','LineWidth',3,'Color',[0 1 0]);
annotation('ellipse',[0.195375 0.288265306122449 0.104 0.211734693877547],'LineStyle','--','LineWidth',3,'Color',[0 1 0]);
annotation('ellipse',[0.347875000000001 0.233826530612243 0.104 0.211734693877547],'LineStyle','--','LineWidth',3,'Color',[0 1 0]);
annotation('ellipse',[0.140625 0.121173469387755 0.0562500000000001 0.0982142857142857],'LineStyle','--','LineWidth',3,'Color',[0 1 0]);
annotation('ellipse',[0.12725 0.303571428571428 0.2115 0.443877551020406],'LineStyle','--','LineWidth',3,'Color',[0 1 1]);
annotation('ellipse',[0.340625 0.21892857142857 0.117500000000001 0.242806122448981],'LineStyle','--','LineWidth',3,'Color',[0 1 1]);
annotation('ellipse',[0.1325 0.107551020408163 0.07 0.123316326530612],'LineStyle','--','LineWidth',3,'Color',[0 1 1]);
hold off

subplot(5,2,4:2:10)
clustering = dendrogram(Tree);
set(clustering,'LineWidth',5)
xlabel('observation index','Fontsize',14),ylabel('joining distance','Fontsize',14)
title({'Dendrogram   (based on distance matrix, D)'; 'and 2 possible horizontal cuts giving 3 or 4 clusters'},'Fontsize',14)
hold on
annotation('line',[0.570625 0.906875],[0.399510204081633 0.400510204081633],'LineStyle','--','LineWidth',3,'Color',[0 1 0]);
annotation('line',[0.570625 0.90625],[0.597214285714286 0.598214285714286],'LineStyle','--','LineWidth',3,'Color',[0 1 1]);
hold off

%------------------------------------------------------------------------------------------------------------------------------------------------------

figure
subplot(5,4,1:4)
axis off
text(-0.05,1,'Comparing inter-cluster distance metrics (linkage methods) for hierarchical agglomerative clustering example','Fontsize',20)
text(0.53,0.75,'-  using Euclidean point-to-point distance metric','Fontsize',18)
text(0.7,0.1,'Dendrograms:','Fontsize',18)

subplot(5,4,[5:6 9:10 13:14 17:18])
plot(X(:,1),X(:,2),'ro','LineWidth',5),xlim([0 10]),ylim([0 10]),grid,axis square,xlabel('x_1','Fontsize',14),ylabel('x_2','Fontsize',14)
title({'10 bivariate observations and'; '2 possible clusterings based on dendrogram cuts'},'Fontsize',14)
hold on
for j = 1:n
    text(X(j,1)+0.2,X(j,2),num2str(j),'Fontsize',14)
end
annotation('ellipse',[0.155375 0.537397959183669 0.104 0.211734693877547],'LineStyle','--','LineWidth',3,'Color',[0 1 0]);
annotation('ellipse',[0.195375 0.288265306122449 0.104 0.211734693877547],'LineStyle','--','LineWidth',3,'Color',[0 1 0]);
annotation('ellipse',[0.347875000000001 0.233826530612243 0.104 0.211734693877547],'LineStyle','--','LineWidth',3,'Color',[0 1 0]);
annotation('ellipse',[0.140625 0.121173469387755 0.0562500000000001 0.0982142857142857],'LineStyle','--','LineWidth',3,'Color',[0 1 0]);
annotation('ellipse',[0.12725 0.303571428571428 0.2115 0.443877551020406],'LineStyle','--','LineWidth',3,'Color',[0 1 1]);
annotation('ellipse',[0.340625 0.21892857142857 0.117500000000001 0.242806122448981],'LineStyle','--','LineWidth',3,'Color',[0 1 1]);
annotation('ellipse',[0.1325 0.107551020408163 0.07 0.123316326530612],'LineStyle','--','LineWidth',3,'Color',[0 1 1]);
annotation('ellipse',[0.344375 0.224489795918367 0.109375 0.230867346938774],'LineStyle','--','LineWidth',3,'Color',[1 0 1]);
annotation('ellipse',[0.1285 0.102040816326531 0.19025 0.404744897959179],'LineStyle','--','LineWidth',3,'Color',[1 0 1]);
annotation('ellipse',[0.146625 0.525510204081631 0.1215 0.238520408163265],'LineStyle','--','LineWidth',3,'Color',[1 0 1]);
hold off

subplot(5,4,[7 11])
Tree = linkage(D,'single');
clustering = dendrogram(Tree);
set(clustering,'LineWidth',5)
ylabel('joining distance','Fontsize',12)
title('single linkage:','Fontsize',14)
hold on
annotation('line',[0.543125 0.699375],[0.589561224489796 0.590561224489796],'LineStyle','--','LineWidth',3,'Color',[0 1 0]);
annotation('line',[0.543125 0.699375],[0.679255102040814 0.680255102040814],'LineStyle','--','LineWidth',3,'Color',[0 1 1]);
hold off

subplot(5,4,[8 12])
Tree = linkage(D,'complete');
clustering = dendrogram(Tree);
set(clustering,'LineWidth',5)
title('complete linkage:','Fontsize',14)
hold on
annotation('line',[0.749375 0.905625],[0.572520408163264 0.573520408163264],'LineStyle','--','LineWidth',3,'Color',[0 1 0]);
annotation('line',[0.748750000000001 0.905000000000001],[0.66563265306122 0.66663265306122],'LineStyle','--','LineWidth',3,'Color',[1 0 1]);
hold off

subplot(5,4,[15 19])
Tree = linkage(D,'average');
clustering = dendrogram(Tree);
set(clustering,'LineWidth',5)
xlabel('observation index','Fontsize',12),ylabel('joining distance','Fontsize',12)
title('average linkage:','Fontsize',14)
hold on
annotation('line',[0.543750000000001 0.700000000000001],[0.228540816326528 0.229540816326528],'LineStyle','--','LineWidth',3,'Color',[0 1 0]);
annotation('line',[0.543125000000001 0.699375000000001],[0.330173469387751 0.331173469387751],'LineStyle','--','LineWidth',3,'Color',[0 1 1]);
hold off

subplot(5,4,[16 20])
Tree = linkage(D,'centroid');
clustering = dendrogram(Tree);
set(clustering,'LineWidth',5)
xlabel('observation index','Fontsize',12)
title('centroid linkage:','Fontsize',14)
hold on
annotation('line',[0.749375000000001 0.905625000000001],[0.232826530612242 0.233826530612242],'LineStyle','--','LineWidth',3,'Color',[0 1 0]);
annotation('line',[0.748750000000001 0.905000000000001],[0.329306122448974 0.330306122448974],'LineStyle','--','LineWidth',3,'Color',[0 1 1]);
hold off

%---------------------------------------------------------------------------------------------------------------------------------
