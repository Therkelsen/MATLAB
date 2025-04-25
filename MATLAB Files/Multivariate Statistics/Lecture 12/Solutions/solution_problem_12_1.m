clear;close all;clc;
format compact

D = [4 6 1 6 9 7 3 10 5 8];

figure
subplot(5,4,1:4)
axis off
text(-0.05,1,'Comparing inter-cluster distance metrics (linkage methods) for hierarchical agglomerative clustering example','Fontsize',20)
text(0.53,0.75,'-  using Euclidean point-to-point distance metric','Fontsize',18)
text(0.7,0.1,'Dendrograms:','Fontsize',18)

subplot(5,4,[7 11])
Tree = linkage(D,'single');
clustering = dendrogram(Tree);
set(clustering,'LineWidth',5)
ylabel('joining distance','Fontsize',12)
title('single linkage:','Fontsize',14)

subplot(5,4,[8 12])
Tree = linkage(D,'complete');
clustering = dendrogram(Tree);
set(clustering,'LineWidth',5)
title('complete linkage:','Fontsize',14)

subplot(5,4,[15 19])
Tree = linkage(D,'average');
clustering = dendrogram(Tree);
set(clustering,'LineWidth',5)
xlabel('observation index','Fontsize',12),ylabel('joining distance','Fontsize',12)
title('average linkage:','Fontsize',14)

subplot(5,4,[16 20])
Tree = linkage(D,'centroid');
clustering = dendrogram(Tree);
set(clustering,'LineWidth',5)
xlabel('observation index','Fontsize',12)
title('centroid linkage:','Fontsize',14)

%---------------------------------------------------------------------------------------------------------------------------------
