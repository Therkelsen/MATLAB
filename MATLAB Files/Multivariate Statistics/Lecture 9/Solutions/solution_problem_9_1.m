clc,clear,close all
format compact
%----------------------------------------
% data & visualization
%----------------------------------------
load('C:\Users\cv\OneDrive - Syddansk Universitet\M-drive moved by SDU-IT\Statistik\MULTISTAT\Lektion 09  Non-parametric statistics\Course material\dataset_problem_9_1.mat')
disp('-----------------------------------')
disp('sign test')
disp('-----------------------------------')
n = length(X)
figure(1)
subplot(121)
edges = 0:10;
histogram(X,edges),title('histogram')
subplot(122)
qqplot(X)
%----------------------------------------
% H0: m = m0
%----------------------------------------
disp('-----------------------------------')
disp('data median and H0 median')
disp('-----------------------------------')
X_median = median(X)
m0 = 4
disp('-----------------------------------')
%----------------------------------------
% test statistic S = (#X > m0) ~ Bin(n,0.5)
%----------------------------------------
disp('test statistic S')
disp('-----------------------------------')
S = length(find(X > m0))
disp('-----------------------------------')
%----------------------------------------
% exact test (binomial distribution)
%----------------------------------------
disp('excact test using binomial')
disp('-----------------------------------')
p_value = 2*min([binocdf(S,n,0.5) 1-binocdf(S,n,0.5)])
%----------------------------------------
% exact test with built-in command
%----------------------------------------
disp('-----------------------------------')
disp('exact test using built-in command')
disp('-----------------------------------')
[p,h,stats] = signtest(X,m0,'method','exact')
%----------------------------------------
% approximate test (standard normal)
%----------------------------------------
disp('----------------------------------------')
disp('approximate test using standard normal')
disp('----------------------------------------')
S_mean = n*0.5;
S_var = n*0.5*(1-0.5);
npos = S;
nneg = n-S;
continuity_correction = 0.5*sign(npos - nneg);
S_Z0 = (S - S_mean - continuity_correction)/sqrt(S_var)
p_value = 2*(1 - normcdf(abs(S_Z0)))
%----------------------------------------
% approximate test with built-in command
%----------------------------------------
disp('----------------------------------------')
disp('approximate test using built-in command')
disp('----------------------------------------')
[p,h,stats] = signtest(X,m0,'method','approximate')
disp('-----------------------------------')
%----------------------------------------

disp('    ')
disp('-----------------------------------')
disp('Wilcoxon signed rank test')
disp('-----------------------------------')
n
disp('-----------------------------------')
%----------------------------------------
% H0: m = m0
%----------------------------------------
disp('data median and H0 median')
disp('---------------------------------------------')
X_median = median(X)
m0 = 4
disp('---------------------------------------------')
%----------------------------------------
% test statistic W
%----------------------------------------
disp('test statistic W')
disp('---------------------------------------------')
abs_dev = abs(X - m0);
signs = sign(X - m0);
[sorted_abs_dev, index] = sort(abs_dev);
W = 0;
for ranks = 1:n
    W = W + ranks*signs(index(ranks));
end
W
%----------------------------------------
% approximate test with W (standard normal)
%----------------------------------------
disp('---------------------------------------------')
disp('approximate test using W and standard normal')
disp('---------------------------------------------')
W_mean = 0;
W_var = n*(n+1)*(2*n+1)/6;
W_Z0 = (W - W_mean)/sqrt(W_var)
p_value = 2*(1 - normcdf(abs(W_Z0)))
%----------------------------------------
% alternative test statistic T (MATLAB)
% (using only smallest of positive
%  and negative ranksums)
% giving same p-value as W
%----------------------------------------
% test statistic T 
%----------------------------------------
disp('---------------------------------------------')
disp('test statistic T (smallest pos/neg ranksum)')
disp('---------------------------------------------')
T_pos = 0;
T_neg = 0;
for ranks = 1:n
    if signs(index(ranks)) > 0
        T_pos = T_pos + ranks;
    else
        T_neg = T_neg + ranks;
    end
end
T = min([T_pos T_neg])
%----------------------------------------
% approximate test with T (standard normal)
%----------------------------------------
disp('---------------------------------------------')
disp('approximate test using T and standard normal')
disp('---------------------------------------------')
T_mean = n*(n+1)/4;
T_var = n*(n+1)*(2*n+1)/24;
T_Z0 = (T - T_mean)/sqrt(T_var)
p_value = 2*(1 - normcdf(abs(T_Z0)))
%----------------------------------------
% approximate test with built-in command
%----------------------------------------
disp('---------------------------------------------')
disp('approximate test with built-in command')
disp('---------------------------------------------')
[p,h,stats] = signrank(X,m0,'method','approximate')
disp('---------------------------------------------')
%----------------------------------------