clear;close all;clc;
format compact

load('M:\Statistik\MULTISTAT\Lektion 12  Cluster Analysis\Course material\dataset_problem_12_3')
load('M:\Statistik\MULTISTAT\Lektion 12  Cluster Analysis\Course material\university_names_problem_12_3')
Z = zscore(X);

disp('-------------------------------------------------------------------------------------------------------')
[n p] = size(Z)
disp('-------------------------------------------------------------------------------------------------------')
D = squareform(pdist(Z,'euclidean'))
disp('-------------------------------------------------------------------------------------------------------')

figure
subplot(5,4,1:4)
axis off
text(-0.05,1,'Comparing inter-cluster distance metrics (linkage methods) for hierarchical agglomerative clustering example','Fontsize',20)
text(0.53,0.75,'-  using Euclidean point-to-point distance metric','Fontsize',18)
text(0.42,0.1,'Raw dendrograms:','Fontsize',18)

subplot(5,4,[5 6 9 10])
Tree = linkage(D,'single');
clustering = dendrogram(Tree);
set(clustering,'LineWidth',5)
ylabel('joining distance','Fontsize',12)
title('single linkage:','Fontsize',14)

subplot(5,4,[7 8 11 12])
Tree = linkage(D,'complete');
clustering = dendrogram(Tree);
set(clustering,'LineWidth',5)
title('complete linkage:','Fontsize',14)

subplot(5,4,[13 14 17 18])
Tree = linkage(D,'average');
clustering = dendrogram(Tree);
set(clustering,'LineWidth',5)
xlabel('observation index','Fontsize',12),ylabel('joining distance','Fontsize',12)
title('average linkage:','Fontsize',14)

subplot(5,4,[15 16 19 20])
Tree = linkage(D,'centroid');
clustering = dendrogram(Tree);
set(clustering,'LineWidth',5)
xlabel('observation index','Fontsize',12)
title('centroid linkage:','Fontsize',14)

%----------------------------------------------------------------------------------------

figure
subplot(5,4,1:4)
axis off
text(-0.05,1,'Comparing inter-cluster distance metrics (linkage methods) for hierarchical agglomerative clustering example','Fontsize',20)
text(0.53,0.75,'-  using Euclidean point-to-point distance metric','Fontsize',18)
text(0.3,0.1,'Possible clusterings using dendrograms:','Fontsize',18)

subplot(5,4,[5 6 9 10])
Tree = linkage(D,'single');
disp('  ')
disp('-------------------------------------------------------------------------------------------------------')
disp('assigments with K = 4 clusters for single linkage:')
disp('-------------------------------------------------------------------------------------------------------')
K = 4;
assignments_clusters_single_linkage = cluster(Tree,'maxclust',K);
cluster_1_id = find(assignments_clusters_single_linkage == 1)'
cluster_1_names = names(cluster_1_id(1:end))
disp('-----------------------------------')
cluster_2_id = find(assignments_clusters_single_linkage == 2)'
cluster_2_names = names(cluster_2_id(1:end))
disp('-----------------------------------')
cluster_3_id = find(assignments_clusters_single_linkage == 3)'
cluster_3_names = names(cluster_3_id(1:end))
disp('-----------------------------------')
cluster_4_id = find(assignments_clusters_single_linkage == 4)'
cluster_4_names = names(cluster_4_id(1:end))
disp('-----------------------------------')
clustering = dendrogram(Tree);
annotation('line',[0.13 0.493125],[0.705632653061218 0.706632653061218],'LineStyle','--','LineWidth',2,'Color',[1 0 0]);
set(clustering,'LineWidth',5)
ylabel('joining distance','Fontsize',12)
title('single linkage with K = 4:','Fontsize',14)

subplot(5,4,[7 8 11 12])
Tree = linkage(D,'complete');
disp('  ')
disp('-------------------------------------------------------------------------------------------------------')
disp('assigments with K = 6 clusters for complete linkage:')
disp('-------------------------------------------------------------------------------------------------------')
K = 6;
assignments_clusters_complete_linkage = cluster(Tree,'maxclust',K);
cluster_1_id = find(assignments_clusters_complete_linkage == 1)'
cluster_1_names = names(cluster_1_id(1:end))
disp('-----------------------------------')
cluster_2_id = find(assignments_clusters_complete_linkage == 2)'
cluster_2_names = names(cluster_2_id(1:end))
disp('-----------------------------------')
cluster_3_id = find(assignments_clusters_complete_linkage == 3)'
cluster_3_names = names(cluster_3_id(1:end))
disp('-----------------------------------')
cluster_4_id = find(assignments_clusters_complete_linkage == 4)'
cluster_4_names = names(cluster_4_id(1:end))
disp('-----------------------------------')
cluster_5_id = find(assignments_clusters_complete_linkage == 5)'
cluster_5_names = names(cluster_5_id(1:end))
disp('-----------------------------------')
cluster_6_id = find(assignments_clusters_complete_linkage == 6)'
cluster_6_names = names(cluster_6_id(1:end))
disp('-----------------------------------')
clustering = dendrogram(Tree);
annotation('line',[0.543125000000001 0.906250000000001],[0.545326530612235 0.546326530612235],'LineStyle','--','LineWidth',2,'Color',[1 0 0]);
set(clustering,'LineWidth',5)
title('complete linkage with K = 6:','Fontsize',14)

subplot(5,4,[13 14 17 18])
Tree = linkage(D,'average');
disp('  ')
disp('-------------------------------------------------------------------------------------------------------')
disp('assigments with K = 5 clusters for average linkage:')
disp('-------------------------------------------------------------------------------------------------------')
K = 5;
assignments_clusters_average_linkage = cluster(Tree,'maxclust',K);
cluster_1_id = find(assignments_clusters_average_linkage == 1)'
cluster_1_names = names(cluster_1_id(1:end))
disp('-----------------------------------')
cluster_2_id = find(assignments_clusters_average_linkage == 2)'
cluster_2_names = names(cluster_2_id(1:end))
disp('-----------------------------------')
cluster_3_id = find(assignments_clusters_average_linkage == 3)'
cluster_3_names = names(cluster_3_id(1:end))
disp('-----------------------------------')
cluster_4_id = find(assignments_clusters_average_linkage == 4)'
cluster_4_names = names(cluster_4_id(1:end))
disp('-----------------------------------')
cluster_5_id = find(assignments_clusters_average_linkage == 5)'
cluster_5_names = names(cluster_5_id(1:end))
disp('-----------------------------------')
clustering = dendrogram(Tree);
annotation('line',[0.131250000000001 0.494375000000001],[0.21920408163264 0.22020408163264],'LineStyle','--','LineWidth',2,'Color',[1 0 0]);
set(clustering,'LineWidth',5)
xlabel('observation index','Fontsize',12),ylabel('joining distance','Fontsize',12)
title('average linkage with K = 5:','Fontsize',14)

subplot(5,4,[15 16 19 20])
Tree = linkage(D,'centroid');
disp('  ')
disp('-------------------------------------------------------------------------------------------------------')
disp('assigments with K = 5 clusters for centroid linkage:')
disp('-------------------------------------------------------------------------------------------------------')
K = 5;
assignments_clusters_centroid_linkage = cluster(Tree,'maxclust',K);
cluster_1_id = find(assignments_clusters_centroid_linkage == 1)'
cluster_1_names = names(cluster_1_id(1:end))
disp('-----------------------------------')
cluster_2_id = find(assignments_clusters_centroid_linkage == 2)'
cluster_2_names = names(cluster_2_id(1:end))
disp('-----------------------------------')
cluster_3_id = find(assignments_clusters_centroid_linkage == 3)'
cluster_3_names = names(cluster_3_id(1:end))
disp('-----------------------------------')
cluster_4_id = find(assignments_clusters_centroid_linkage == 4)'
cluster_4_names = names(cluster_4_id(1:end))
disp('-----------------------------------')
cluster_5_id = find(assignments_clusters_centroid_linkage == 5)'
cluster_5_names = names(cluster_5_id(1:end))
disp('-----------------------------------')
clustering = dendrogram(Tree);
annotation('line',[0.543750000000001 0.906875000000001],[0.214510204081619 0.215510204081619],'LineStyle','--','LineWidth',2,'Color',[1 0 0]);
set(clustering,'LineWidth',5)
xlabel('observation index','Fontsize',12)
title('centroid linkage with K = 5:','Fontsize',14)

%----------------------------------------------------------------------------------------

figure
subplot(6,4,1:4)
axis off
text(0,1.2,'Non-hierarchical "K-means" clustering example:','Fontsize',24)
text(0,0.85,'n=25 observations, p=6 variables, K=3,4,5,6 clusters','Fontsize',18)
text(0,0.55,'SSW minimized over 10 initial clusterings pr. K-value','Fontsize',16)
text(0.65,1.2,'Silhouette value pr. obs  j = 1,...,n','Fontsize',12)
text(0.65,0.9,'h_s_,_j = (d_a_v_e_r_,_n_e_x_t_ _c_l_u_s_t_e_r - d_a_v_e_r_,_o_w_n_ _c_l_u_s_t_e_r) / max(these 2 numbers)','Fontsize',12)
text(0.65,0.55,'=>  -1 < h_s_,_j < 1','Fontsize',12)

disp('  ')
disp('-------------------------------------------------------------------------------------------------------')
disp('assigments with K = 3 clusters for Kmeans algorithm:')
disp('-------------------------------------------------------------------------------------------------------')
K = 3;
[assignments_clusters_Kmeans_3, C, cluster_SSWs] = kmeans(Z,K,'distance','sqEuclidean','start','sample','replicates',10);
cluster_SSWs
total_cluster_SSW = sum(cluster_SSWs)
cluster_1_id = (find(assignments_clusters_Kmeans_3 == 1))'
cluster_1_names = names(cluster_1_id(1:end))
cluster_2_id = (find(assignments_clusters_Kmeans_3 == 2))'
cluster_2_names = names(cluster_2_id(1:end))
cluster_3_id = (find(assignments_clusters_Kmeans_3 == 3))'
cluster_3_names = names(cluster_3_id(1:end))
disp('-----------------------------------')
subplot(6,4,[6 10])
[silh h] = silhouette(Z,assignments_clusters_Kmeans_3);title(['mean(silhouette) = ' num2str(mean(silh),3)],'Fontsize',14)

disp('  ')
disp('-------------------------------------------------------------------------------------------------------')
disp('assigments with K = 4 clusters for Kmeans algorithm:')
disp('-------------------------------------------------------------------------------------------------------')
K = 4;
[assignments_clusters_Kmeans_4, C, cluster_SSWs] = kmeans(Z,K,'distance','sqEuclidean','start','sample','replicates',10);
cluster_SSWs
total_cluster_SSW = sum(cluster_SSWs)
cluster_1_id = (find(assignments_clusters_Kmeans_4 == 1))'
cluster_1_names = names(cluster_1_id(1:end))
cluster_2_id = (find(assignments_clusters_Kmeans_4 == 2))'
cluster_2_names = names(cluster_2_id(1:end))
cluster_3_id = (find(assignments_clusters_Kmeans_4 == 3))'
cluster_3_names = names(cluster_3_id(1:end))
cluster_4_id = (find(assignments_clusters_Kmeans_4 == 4))'
cluster_4_names = names(cluster_4_id(1:end))
disp('-----------------------------------')
subplot(6,4,[8 12])
[silh h] = silhouette(Z,assignments_clusters_Kmeans_4);title(['mean(silhouette) = ' num2str(mean(silh),3)],'Fontsize',14)

disp('  ')
disp('-------------------------------------------------------------------------------------------------------')
disp('assigments with K = 5 clusters for Kmeans algorithm:')
disp('-------------------------------------------------------------------------------------------------------')
K = 5;
[assignments_clusters_Kmeans_5, C, cluster_SSWs] = kmeans(Z,K,'distance','sqEuclidean','start','sample','replicates',10);
cluster_SSWs
total_cluster_SSW = sum(cluster_SSWs)
cluster_1_id = (find(assignments_clusters_Kmeans_5 == 1))'
cluster_1_names = names(cluster_1_id(1:end))
cluster_2_id = (find(assignments_clusters_Kmeans_5 == 2))'
cluster_2_names = names(cluster_2_id(1:end))
cluster_3_id = (find(assignments_clusters_Kmeans_5 == 3))'
cluster_3_names = names(cluster_3_id(1:end))
cluster_4_id = (find(assignments_clusters_Kmeans_5 == 4))'
cluster_4_names = names(cluster_4_id(1:end))
cluster_5_id = (find(assignments_clusters_Kmeans_5 == 5))'
cluster_5_names = names(cluster_5_id(1:end))
disp('-----------------------------------')
subplot(6,4,[18 22])
[silh h] = silhouette(Z,assignments_clusters_Kmeans_5);title(['mean(silhouette) = ' num2str(mean(silh),3)],'Fontsize',14)

disp('  ')
disp('-------------------------------------------------------------------------------------------------------')
disp('assigments with K = 6 clusters for Kmeans algorithm:')
disp('-------------------------------------------------------------------------------------------------------')
K = 6;
[assignments_clusters_Kmeans_6, C, cluster_SSWs] = kmeans(Z,K,'distance','sqEuclidean','start','sample','replicates',10);
cluster_SSWs
total_cluster_SSW = sum(cluster_SSWs)
cluster_1_id = (find(assignments_clusters_Kmeans_6 == 1))'
cluster_1_names = names(cluster_1_id(1:end))
cluster_2_id = (find(assignments_clusters_Kmeans_6 == 2))'
cluster_2_names = names(cluster_2_id(1:end))
cluster_3_id = (find(assignments_clusters_Kmeans_6 == 3))'
cluster_3_names = names(cluster_3_id(1:end))
cluster_4_id = (find(assignments_clusters_Kmeans_6 == 4))'
cluster_4_names = names(cluster_4_id(1:end))
cluster_5_id = (find(assignments_clusters_Kmeans_6 == 5))'
cluster_5_names = names(cluster_5_id(1:end))
cluster_6_id = (find(assignments_clusters_Kmeans_6 == 6))'
cluster_6_names = names(cluster_6_id(1:end))
disp('-----------------------------------')
subplot(6,4,[20 24])
[silh h] = silhouette(Z,assignments_clusters_Kmeans_6);title(['mean(silhouette) = ' num2str(mean(silh),3)],'Fontsize',14)
%---------------------------------------------------------------------------------------------------------------------------------
