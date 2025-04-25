clc,clear,close all
format compact
%----------------------------------------
% data 
%----------------------------------------
load('C:\Users\cv\OneDrive - Syddansk Universitet\M-drive moved by SDU-IT\Statistik\MULTISTAT\Lektion 09  Non-parametric statistics\Course material\dataset_problem_9_2.mat')
disp('-----------------------------------')
disp('Mann-Whitney U test')
disp('-----------------------------------')
nx = length(X)
ny = length(Y)
%----------------------------------------
% visualization of data
%----------------------------------------
figure(1)
subplot(221)
edges = 0:1:20;
histogram(X,edges),title('histogram X'),xlim([0 20])
subplot(222)
qqplot(X),title('QQplot X')
subplot(223)
histogram(Y,edges),title('histogram Y'),xlim([0 20])
subplot(224)
qqplot(Y),title('QQplot Y')
%----------------------------------------
% H0: mx = my
%----------------------------------------
disp('---------------------------------------------')
disp('data medians')
disp('---------------------------------------------')
X_median = median(X)
Y_median = median(Y)
disp('---------------------------------------------')
%----------------------------------------
% test statistic U
%----------------------------------------
disp('test statistic U (Mann-Whitney)')
disp('---------------------------------------------')
U = 0;
for i = 1:nx
    for j = 1:ny
        if X(i) > Y(j)
            U = U + 1;
        else if X(i) == Y(j)
                U = U + 1/2;
            end
        end
    end
end
U
%----------------------------------------
% approximate test with U (standard normal)
%----------------------------------------
disp('---------------------------------------------')
disp('approximate test using U and standard normal')
disp('---------------------------------------------')
U_mean = nx*ny/2;
U_var = nx*ny*(nx+ny+1)/12;
continuity_correction = 0.5*sign(U - U_mean);
U_Z0 = (U - U_mean - continuity_correction)/sqrt(U_var)
p_value = 2*(1 - normcdf(abs(U_Z0)))
%----------------------------------------
% alternative test statistic W (MATLAB)
% giving same p-value as U
% (W is the ranksum for X, i.e. the sum of
%  all ranks for X when X and Y data are
%  are ranked together from 1 to nx+ny)
%----------------------------------------
% test statistic W
%----------------------------------------
disp('---------------------------------------------')
disp('test statistic W (ranksum for X data)')
disp('---------------------------------------------')
W = U + nx*(nx+1)/2
%----------------------------------------
% approximate test with W (standard normal)
%----------------------------------------
disp('---------------------------------------------')
disp('approximate test using W and standard normal')
disp('---------------------------------------------')
W_mean = U_mean + nx*(nx+1)/2;
W_var = U_var;
continuity_correction = 0.5*sign(W - W_mean);
W_Z0 = (W - W_mean - continuity_correction)/sqrt(W_var)
p_value = 2*(1 - normcdf(abs(W_Z0)))
%----------------------------------------
% approximate test with built-in command
%----------------------------------------
disp('---------------------------------------------')
disp('approximate test with built-in command')
disp('---------------------------------------------')
[p,h,stats] = ranksum(X,Y,'method','approximate')
disp('---------------------------------------------')
%----------------------------------------