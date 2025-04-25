clc,clear,close all
format compact
%----------------------------------------
% data 
%----------------------------------------
load('C:\Users\cv\OneDrive - Syddansk Universitet\M-drive moved by SDU-IT\Statistik\MULTISTAT\Lektion 09  Non-parametric statistics\Course material\dataset_problem_9_3a.mat')
disp('-----------------------------------')
disp('Kruskal-Wallis test')
disp('-----------------------------------')
n1 = length(X1);
n2 = length(X2);
n3 = length(X3);
%----------------------------------------
g = 3
n1, n2, n3
N = n1 + n2 + n3
%----------------------------------------
% visualization of data
%----------------------------------------
figure(1)
subplot(321)
edges = 0:1:20;
histogram(X1,edges),title('histogram X1'),xlim([0 20])
subplot(322)
qqplot(X1),title('QQplot X1')
subplot(323)
histogram(X2,edges),title('histogram X2'),xlim([0 20])
subplot(324)
qqplot(X2),title('QQplot X2')
subplot(325)
histogram(X3,edges),title('histogram X3'),xlim([0 20])
subplot(326)
qqplot(X3),title('QQplot X3')
%----------------------------------------
% H0: m1 = m2 = m3
%----------------------------------------
disp('---------------------------------------------')
disp('data medians')
disp('---------------------------------------------')
X1_median = median(X1)
X2_median = median(X2)
X3_median = median(X3)
disp('---------------------------------------------')
%----------------------------------------
% test statistic H
%----------------------------------------
disp('test statistic H (Kruskall-Wallis)')
disp('---------------------------------------------')
X = [X1 X2 X3];
group = [1*ones(1,n1) 2*ones(1,n2) 3*ones(1,n3)];
[X_sorted, index] = sort(X);
group_for_X_sorted = group(index);
ranks = 1:N;
mean_rank = mean(ranks)
%----------------------------------------
n = [n1 n2 n3];
Hnum = 0;
for i = 1:g
    group_mean_rank = mean(ranks(find(group_for_X_sorted == i)))
    Hnum = Hnum + n(i)*(group_mean_rank - mean_rank)^2;
end
Hnum = (N-1)*Hnum;
SSB = Hnum
%----------------------------------------
Hden = 0;
for i = 1:g
    group_ranks = ranks(find(group_for_X_sorted == i));
    for j = 1:n(i)
        Hden = Hden + (group_ranks(j) - mean_rank)^2;
    end
end
SST = Hden
%----------------------------------------
H = SSB/SST
%----------------------------------------
% approximate test with H (chi-squared)
%----------------------------------------
disp('-----------------------------------------------------')
disp('approximate test using H and chi-squared distribution')
disp('-----------------------------------------------------')
df = g-1
p_value = 1-chi2cdf(H,df)
%----------------------------------------
% approximate test with built-in command
%----------------------------------------
disp('---------------------------------------------')
disp('approximate test with built-in command')
disp('---------------------------------------------')
[p,tbl,stats] = kruskalwallis(X,group)
disp('---------------------------------------------')

%------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------

disp('    ')
%----------------------------------------
% data 
%----------------------------------------
load('C:\Users\cv\OneDrive - Syddansk Universitet\M-drive moved by SDU-IT\Statistik\MULTISTAT\Lektion 09  Non-parametric statistics\Course material\dataset_problem_9_3b.mat')
disp('-----------------------------------')
disp('Friedman test')
disp('-----------------------------------')
%----------------------------------------
g = 4
n = 50
N = g*n
%----------------------------------------
% visualization of data
%----------------------------------------
figure(4)
subplot(221)
edges = 0:1:20;
histogram(X1,edges),title('histogram X1'),xlim([0 20])
subplot(222)
histogram(X2,edges),title('histogram X2'),xlim([0 20])
subplot(223)
histogram(X3,edges),title('histogram X3'),xlim([0 20])
subplot(224)
histogram(X4,edges),title('histogram X4'),xlim([0 20])
%----------------------------------------
figure(5)
subplot(221)
qqplot(X1),title('QQplot X1')
subplot(222)
qqplot(X2),title('QQplot X2')
subplot(223)
qqplot(X3),title('QQplot X3')
subplot(224)
qqplot(X4),title('QQplot X4')
%----------------------------------------
% H0: m1 = m2 = m3 = m4
%----------------------------------------
disp('---------------------------------------------')
disp('data medians')
disp('---------------------------------------------')
X1_median = median(X(:,1))
X2_median = median(X(:,2))
X3_median = median(X(:,3))
X4_median = median(X(:,4))
disp('---------------------------------------------')
%----------------------------------------
% test statistic Q
%----------------------------------------
disp('test statistic Q (Friedman)')
disp('---------------------------------------------')
X = [X1 X2 X3 X4];
ranks_matrix = zeros(n,g);
for j = 1:n
    [row_rank,TIEADJ] = tiedrank(X(j,:));
    ranks_matrix(j,:) = row_rank;
end
column_mean_ranks = mean(ranks_matrix)
Q = 0;
for i = 1:g
    Q = Q + (12*n/(g*(g+1)))*(column_mean_ranks(i) - (g+1)/2)^2;
end
Q
%----------------------------------------
% approximate test with Q (chi-squared)
%----------------------------------------
disp('-----------------------------------------------------')
disp('approximate test using Q and chi-squared distribution')
disp('-----------------------------------------------------')
df = g-1
p_value = 1-chi2cdf(Q,df)
%----------------------------------------
% approximate test with built-in command
%----------------------------------------
disp('---------------------------------------------')
disp('approximate test with built-in command')
disp('---------------------------------------------')
reps = 1
[p,tbl,stats] = friedman(X,reps)
disp('---------------------------------------------')
%----------------------------------------
