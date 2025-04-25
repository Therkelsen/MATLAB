clc,clear,close all
format compact
%----------------------------------------
% data 
%----------------------------------------
load('C:\Users\cv\OneDrive - Syddansk Universitet\M-drive moved by SDU-IT\Statistik\MULTISTAT\Lektion 09  Non-parametric statistics\Course material\dataset_problem_9_6.mat')
disp('-----------------------------------')
disp('Contingency table test')
disp('-----------------------------------')
X_obs
rows = size(X_obs,1)
cols = size(X_obs,2)
N = sum(sum(X_obs))
%----------------------------------------
% test statistic T (Pearson's X^2)
%----------------------------------------
disp('---------------------------------------------')
disp('test statistic T (Pearson''s X^2)')
disp('---------------------------------------------')
row_sums = sum(X_obs,2)
column_sums = sum(X_obs,1)
X_exp = zeros(rows,cols);
for i = 1:rows
    for j = 1:cols
        X_exp(i,j) = sum(X_obs(i,:))*sum(X_obs(:,j))/sum(sum(X_obs));
    end
end
X_exp
T = sum(sum((X_obs - X_exp).^2 ./ X_exp))
%----------------------------------------
% approximate test with T (Pearson''s GOF X^2)
%----------------------------------------
disp('-----------------------------------------------------')
disp('approximate test with T (Pearson''s GOF X^2)')
disp('-----------------------------------------------------')
alpha = 0.05
df = (rows-1)*(cols-1)
critical_value = chi2inv(1-alpha,df)
reject_H0 = T > critical_value
p_value = 1-chi2cdf(T,df)
%----------------------------------------
% approximate test with built-in command
%----------------------------------------
disp('---------------------------------------------')
disp('approximate test with built-in command (Pearson''s GOF X^2)')
disp('---------------------------------------------')
Y = [];
for i = 1:rows
    for j = 1:cols
        Y = [Y; repmat([i j],X_obs(i,j),1)];
    end
end
[contingency_table,T,p_value] = crosstab(Y(:,1),Y(:,2))
disp('---------------------------------------------')
%----------------------------------------