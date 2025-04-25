clc,clear,close all
format compact
%----------------------------------------
% data 
%----------------------------------------
load('C:\Users\cv\OneDrive - Syddansk Universitet\M-drive moved by SDU-IT\Statistik\MULTISTAT\Lektion 09  Non-parametric statistics\Course material\dataset_problem_9_4a.mat')
disp('-----------------------------------')
disp('One-sample Kolmogorov-Smirnov test')
disp('-----------------------------------')
n = length(X)
%----------------------------------------
% hypothesized distribution 
%----------------------------------------
disp('---------------------------------------------')
disp('hypothesized distribution')
disp('---------------------------------------------')
lambda_0 = 10;
prob_dist_0 = makedist('Exponential','mu',1/lambda_0)
%----------------------------------------
% visualization of data vs. model
%----------------------------------------
figure(1)
qqplot(X,prob_dist_0),title('QQplot data vs. model'),xlabel('model quantiles'),ylabel('data quantiles')
%------------------------------------
figure(2)
bins = 10;
histogram(X,bins,'Normalization','pdf');
hold on
x = 0:0.01:5/lambda_0;
X_0_pdf = exppdf(x,1/lambda_0);       
plot(x,X_0_pdf,'r','Linewidth',2),title('data vs. model pdf')   
legend(' data histogram',' model distribution pdf','Location','NorthEast')
hold off
%--------------------------------------
figure(3)
H = cdfplot(X);
H.LineWidth = 1.5;
hold on
X_0_cdf = expcdf(x,1/lambda_0);       
plot(x,X_0_cdf,'r','Linewidth',1.5),title('data vs. model cdf')  
%----------------------------------------
% test statistic Dn
%----------------------------------------
disp('---------------------------------------------')
disp('test statistic Dn (Kolmogorov_Smirnov)')
disp('---------------------------------------------')
X_cdf_x_values = H.XData;
X_cdf_F_values = H.YData;
X_0_cdf_F_values = expcdf(X_cdf_x_values,1/lambda_0); 
cdf_diff_values = abs(X_cdf_F_values - X_0_cdf_F_values);
Dn = max(cdf_diff_values)
I = find(cdf_diff_values == Dn);
if X_0_cdf_F_values(I) > X_cdf_F_values(I)
    line([X_cdf_x_values(I)-0.001 X_cdf_x_values(I)-0.001],...
       [min([X_cdf_F_values(I) X_0_cdf_F_values(I)]) max([X_cdf_F_values(I) X_0_cdf_F_values(I)])],...
       'Color','k','LineWidth',2,'LineStyle',':')
else
    line([X_cdf_x_values(I)+0.001 X_cdf_x_values(I)+0.001],...
       [min([X_cdf_F_values(I) X_0_cdf_F_values(I)]) max([X_cdf_F_values(I) X_0_cdf_F_values(I)])],...
       'Color','k','LineWidth',2,'LineStyle',':')
end
legend(' data cdf',' model distribution cdf',' test statistic Dn','Location','East')
hold off
%----------------------------------------
% approximate test with Dn
%----------------------------------------
disp('-----------------------------------------------------')
disp('approximate test')
disp('-----------------------------------------------------')
alpha = 0.05
critical_value = sqrt(-0.5*log(alpha/2))/sqrt(n)
reject_H0 = Dn > critical_value
%----------------------------------------
% approximate test with built-in command
%----------------------------------------
disp('---------------------------------------------')
disp('approximate test with built-in command')
disp('---------------------------------------------')
[reject_H0,p_value,Dn,critical_value] = kstest(X,'CDF',prob_dist_0)
disp('---------------------------------------------')
%----------------------------------------

%------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------

disp('    ')
%----------------------------------------
% data 
%----------------------------------------
load('C:\Users\cv\OneDrive - Syddansk Universitet\M-drive moved by SDU-IT\Statistik\MULTISTAT\Lektion 09  Non-parametric statistics\Course material\dataset_problem_9_4b.mat')
disp('-----------------------------------')
disp('Two-sample Kolmogorov-Smirnov test')
disp('-----------------------------------')
%----------------------------------------
nx = length(X)
ny = length(Y)
%----------------------------------------
% visualization of data X vs data Y
%----------------------------------------
figure(4)
qqplot(X,Y),title('QQplot data X vs. data Y'),xlabel('data Y quantiles'),ylabel('data X quantiles')
%------------------------------------
figure(5)
bins = 10;
subplot(211)
histogram(X,bins,'Normalization','pdf'),title('data X histogram');
subplot(212)
histogram(Y,bins,'Normalization','pdf'),title('data Y histogram');
%--------------------------------------
figure(6)
Hx = cdfplot(X);
Hx.LineWidth = 1.5;
Hx.Color = 'b';
hold on
Hy = cdfplot(Y);
Hy.LineWidth = 1.5;
Hy.Color = 'r';  
%----------------------------------------
% test statistic D
%----------------------------------------
disp('---------------------------------------------')
disp('test statistic D (Kolmogorov_Smirnov)')
disp('---------------------------------------------')
X_cdf_x_values = Hx.XData;
X_cdf_F_values = Hx.YData;
Y_cdf_x_values = Hy.XData;
Y_cdf_F_values = Hy.YData;
xmin = min([X_cdf_x_values(2) Y_cdf_x_values(2)]);
xmax = max([X_cdf_x_values(2*nx+1) Y_cdf_x_values(2*ny+1)]);
dx = (xmax - xmin)/1e3;
x = xmin:dx:xmax;
cdf_abs_diff = zeros(1,length(x));
for i = 1:length(x)
    Ix = find(X_cdf_x_values <= x(i));
    Fx = max(X_cdf_F_values(Ix));
    Iy = find(Y_cdf_x_values <= x(i));
    Fy = max(Y_cdf_F_values(Iy));
    cdf_abs_diff_values(i) = abs(Fx - Fy);
end
D = max(cdf_abs_diff_values)
%-----------------------------------
I = find(cdf_abs_diff_values == D,1);
xline = x(I);
Ix = find(X_cdf_x_values <= xline);
Fx = max(X_cdf_F_values(Ix));
Iy = find(Y_cdf_x_values <= xline);
Fy = max(Y_cdf_F_values(Iy));
if Fx > Fy
    line([xline xline],[Fy Fx],'Color','k','LineWidth',2,'LineStyle',':')
else
    line([xline xline],[Fx Fy],'Color','k','LineWidth',2,'LineStyle',':')
end
legend(' data X cdf',' data Y cdf',' test statistic D','Location','East')
hold off
%----------------------------------------
% approximate test with D
%----------------------------------------
disp('-----------------------------------------------------')
disp('approximate test')
disp('-----------------------------------------------------')
alpha = 0.05
critical_value = sqrt(-0.5*log(alpha/2))*sqrt(1/nx + 1/ny)
reject_H0 = D > critical_value
%----------------------------------------
% approximate test with built-in command
%----------------------------------------
disp('---------------------------------------------')
disp('approximate test with built-in command')
disp('---------------------------------------------')
[reject_H0,p_value,D] = kstest2(X,Y,'Tail','unequal')
disp('---------------------------------------------')
%----------------------------------------

