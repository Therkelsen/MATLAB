clear;close all;clc;
format compact

disp('----------------------------------------------------------------------')
X = [5 4; 1 -2; -1 1; 3 1]
[n p] = size(X)

disp('---------------------------------------------------------------------')
disp('init clusters (12) and (34):')
disp('---------------------------------------------------------------------')
K = 2;
[idx2, cluster_centroids, cluster_SSWs] = kmeans(X,K,'distance','sqEuclidean','start',[3 1; 1 1]);
final_cluster_1 = (find(idx2 == 1))'
final_cluster_2 = (find(idx2 == 2))'
cluster_centroids
cluster_SSWs

disp('---------------------------------------------------------------------')
disp('init clusters (13) and (24):')
disp('---------------------------------------------------------------------')
K = 2;
[idx2, cluster_centroids, cluster_SSWs] = kmeans(X,K,'distance','sqEuclidean','start',[2 2.5; 2 -0.5]);
final_cluster_1 = (find(idx2 == 1))'
final_cluster_2 = (find(idx2 == 2))'
cluster_centroids
cluster_SSWs

disp('---------------------------------------------------------------------')
disp('init clusters (1) and (234):')
disp('---------------------------------------------------------------------')
K = 2;
[idx2, cluster_centroids, cluster_SSWs] = kmeans(X,K,'distance','sqEuclidean','start',[5 4; 1 0]);
final_cluster_1 = (find(idx2 == 1))'
final_cluster_2 = (find(idx2 == 2))'
cluster_centroids
cluster_SSWs

disp('---------------------------------------------------------------------')
disp('init clusters (123) and (4):')
disp('---------------------------------------------------------------------')
K = 2;
[idx2, cluster_centroids, cluster_SSWs] = kmeans(X,K,'distance','sqEuclidean','start',[5/3 1; 3 1]);
final_cluster_1 = (find(idx2 == 1))'
final_cluster_2 = (find(idx2 == 2))'
cluster_centroids
cluster_SSWs
disp('---------------------------------------------------------------------')
%---------------------------------------------------------------------------------------------------------------------------------
