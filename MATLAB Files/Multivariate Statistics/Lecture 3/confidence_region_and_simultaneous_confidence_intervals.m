clc;clear;close all
format compact
%--------------------------------------------------------------------------
% define MVN
%--------------------------------------------------------------------------
mu1 = 4;
mu2 = 3;
sigma1 = 1;
sigma2 = 1;
var1 = sigma1^2;
var2 = sigma2^2;
rho = 0.75;
mu = [mu1 mu2];
Sigma = [var1 rho*sqrt(var1*var2); rho*sqrt(var1*var2) var2];
maxvar = max(var1,var2);
%--------------------------------------------------------------------------
% generate observations from MVN
%--------------------------------------------------------------------------
p = 2;
n = 100;
obs = mvnrnd(mu,Sigma,n);
%--------------------------------------------------------------------------
% descriptive statistics
%--------------------------------------------------------------------------
mu_hat = mean(obs);
mu1_hat = mu_hat(1);
mu2_hat = mu_hat(2);
Sigma_hat = cov(obs);
subplot(4,2,3)
axis off
%--------------------------------------------------------------------------
%  plot observations and prediction ellipses
%--------------------------------------------------------------------------
x1 = mu1-3*sigma1:6*sigma1/50:mu1+3*sigma1; 
x2 = mu2-3*sigma2:6*sigma2/50:mu2+3*sigma2;
[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],mu,Sigma);
F = reshape(F,length(x2),length(x1));
subplot(4,2,[5 7])
alpha = [0.5 0.25 0.05];
plot(obs(:,1),obs(:,2),'+'),grid,axis([mu1-3*sqrt(maxvar) mu1+3*sqrt(maxvar) mu2-3*sqrt(maxvar) mu2+3*sqrt(maxvar)]),...
    xlabel('x_1','Fontsize',16); ylabel('x_2','Fontsize',16),...
    title({['n = ' num2str(n) ' observations scatter plot for (X_1, X_2) and'];['(1-{\alpha})*100% prediction ellipses for  {\alpha}  =  [' num2str(alpha(3)) '  ' num2str(alpha(2)) '  ' num2str(alpha(1)) ']']},'Fontsize',14)
hold on
c = chi2inv(1-alpha,p);
C = (1/(2*pi*sqrt(det(Sigma))))*exp(-c/2);
contour(x1,x2,F,C,'Color','k','LineWidth',2),axis square
text(-4.5,5,'Descriptive statistics:','Fontsize',14)
text(-4.5,4.2,['{\mu}_h_a_t = x = [' num2str(mu1_hat,3) '  ' num2str(mu2_hat,3) ']'],'Fontsize',14)
text(-3.5,4.65,'_','Fontsize',16)
text(-4.5,3.5,['{\Sigma}_h_a_t = S = [' num2str(Sigma_hat(1,1),3) '  ' num2str(Sigma_hat(1,2),3) ';'],'Fontsize',14)
text(-4.5,3.2,['                  ' num2str(Sigma_hat(2,1),3) '  ' num2str(Sigma_hat(2,2),3) ']'],'Fontsize',14)
hold off
%--------------------------------------------------------------------------
% Confidence region for mu (not for observations !)
%--------------------------------------------------------------------------
subplot(4,2,2:2:6)
y1 = mu1_hat-sqrt(Sigma_hat(1,1)):sqrt(Sigma_hat(1,1))/50:mu1_hat+sqrt(Sigma_hat(1,1)); 
y2 = mu2_hat-sqrt(Sigma_hat(2,2)):sqrt(Sigma_hat(2,2))/50:mu2_hat+sqrt(Sigma_hat(2,2));
[Y1,Y2] = meshgrid(y1,y2);
F = mvnpdf([Y1(:) Y2(:)],mu_hat,Sigma_hat);
F = reshape(F,length(y2),length(y1));
alpha = 0.05;
alpha_contour = [0.05 0.05];  % must be a vector for "contour"
c = (p*(n-1)/(n*(n-p)))*finv(1-alpha_contour,p,n-p);
C = (1/(2*pi*sqrt(det(Sigma_hat))))*exp(-c/2);
contour(y1,y2,F,C,'Color','k','LineWidth',3),axis square,grid,xlabel('{\mu}_1','Fontsize',18),ylabel('{\mu}_2','Fontsize',18),...
        title({'95% Confidence region and intervals for ({\mu}_1, {\mu}_2) based on observations:';'CR (black), simult.CI (red), marg.CI (blue), Bonf.CI (green)'},'Fontsize',14)
hold on
[V, lambda] = eig(Sigma_hat);
ellipse_axes_half_lengths_1 = sqrt(lambda(1,1))*sqrt((p*(n-1)/(n*(n-p)))*finv(1-alpha,p,n-p));
ellipse_axes_half_lengths_2 = sqrt(lambda(2,2))*sqrt((p*(n-1)/(n*(n-p)))*finv(1-alpha,p,n-p));
line([mu1_hat  mu1_hat+ellipse_axes_half_lengths_1*V(1,1)],[mu2_hat mu2_hat+ellipse_axes_half_lengths_1*V(2,1)],'LineWidth',2,'Color','k')
line([mu1_hat  mu1_hat+ellipse_axes_half_lengths_2*V(1,2)],[mu2_hat mu2_hat+ellipse_axes_half_lengths_2*V(2,2)],'LineWidth',2,'Color','k')
%--------------------------------------------------------------------------
% Simultaneous confidence intervals for mu (not for observations !)
%--------------------------------------------------------------------------
mu1_sim_CI = [mu1_hat-sqrt((p*(n-1)/(n-p))*finv(1-alpha,p,n-p))*sqrt(Sigma_hat(1,1)/n) mu1_hat+sqrt((p*(n-1)/(n-p))*finv(1-alpha,p,n-p))*sqrt(Sigma_hat(1,1)/n)];
line([mu1_sim_CI(1) mu1_sim_CI(1)],[0 6],'LineWidth',2,'Color','r')
line([mu1_sim_CI(2) mu1_sim_CI(2)],[0 6],'LineWidth',2,'Color','r')
mu2_sim_CI = [mu2_hat-sqrt((p*(n-1)/(n-p))*finv(1-alpha,p,n-p))*sqrt(Sigma_hat(2,2)/n) mu2_hat+sqrt((p*(n-1)/(n-p))*finv(1-alpha,p,n-p))*sqrt(Sigma_hat(2,2)/n)];
line([1 7],[mu2_sim_CI(1) mu2_sim_CI(1)],'LineWidth',2,'Color','r')
line([1 7],[mu2_sim_CI(2) mu2_sim_CI(2)],'LineWidth',2,'Color','r'),xlim(mu1_sim_CI+[-0.1 0.1]),ylim(mu2_sim_CI+[-0.1 0.1])
%--------------------------------------------------------------------------
% Marginal confidence intervals for mu (not for observations !)
%--------------------------------------------------------------------------
mu1_marg_CI = [mu1_hat-tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(1,1)/n) mu1_hat+tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(1,1)/n)];
line([mu1_marg_CI(1) mu1_marg_CI(1)],[0 6],'LineWidth',2,'Color','b')
line([mu1_marg_CI(2) mu1_marg_CI(2)],[0 6],'LineWidth',2,'Color','b')
mu2_marg_CI = [mu2_hat-tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(2,2)/n) mu2_hat+tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(2,2)/n)];
line([1 7],[mu2_marg_CI(1) mu2_marg_CI(1)],'LineWidth',2,'Color','b')
line([1 7],[mu2_marg_CI(2) mu2_marg_CI(2)],'LineWidth',2,'Color','b')
%--------------------------------------------------------------------------
% Bonferroni simultaneous confidence intervals for mu (not for observations !)
%--------------------------------------------------------------------------
mu1_bonf_CI = [mu1_hat-tinv(1-alpha/(2*p),n-1)*sqrt(Sigma_hat(1,1)/n) mu1_hat+tinv(1-alpha/(2*p),n-1)*sqrt(Sigma_hat(1,1)/n)];
line([mu1_bonf_CI(1) mu1_bonf_CI(1)],[0 6],'LineWidth',2,'Color','g')
line([mu1_bonf_CI(2) mu1_bonf_CI(2)],[0 6],'LineWidth',2,'Color','g')
mu2_bonf_CI = [mu2_hat-tinv(1-alpha/(2*p),n-1)*sqrt(Sigma_hat(2,2)/n) mu2_hat+tinv(1-alpha/(2*p),n-1)*sqrt(Sigma_hat(2,2)/n)];
line([1 7],[mu2_bonf_CI(1) mu2_bonf_CI(1)],'LineWidth',2,'Color','g')
line([1 7],[mu2_bonf_CI(2) mu2_bonf_CI(2)],'LineWidth',2,'Color','g')
text(0.999*mu1_hat,mu2_hat,'o','Fontsize',14)
text(mu1_hat-0.075,mu2_hat+0.03,'({\mu}_1_,_h_a_t, {\mu}_2_,_h_a_t)','Fontsize',14)
hold off
%--------------------------------------------------------------------------
% Test confidence region and intervals
%--------------------------------------------------------------------------
N = 1000;
mu_in_CR = 0;
mu_in_sim_CI = 0;
mu_in_marg_CI = 0;
mu_in_bonf_CI = 0;
for iterations = 1:N
    obs = mvnrnd(mu,Sigma,n);
    mu_hat = mean(obs);
    mu1_hat = mu_hat(1);
    mu2_hat = mu_hat(2);
    Sigma_hat = cov(obs);
    if  (mu_hat-mu)*inv(Sigma_hat)*(mu_hat-mu)' <= (p*(n-1)/(n*(n-p)))*finv(1-alpha,p,n-p)
        mu_in_CR = mu_in_CR + 1;
    end
    mu1_sim_CI = [mu1_hat-sqrt((p*(n-1)/(n-p))*finv(1-alpha,p,n-p))*sqrt(Sigma_hat(1,1)/n) mu1_hat+sqrt((p*(n-1)/(n-p))*finv(1-alpha,p,n-p))*sqrt(Sigma_hat(1,1)/n)];
    mu2_sim_CI = [mu2_hat-sqrt((p*(n-1)/(n-p))*finv(1-alpha,p,n-p))*sqrt(Sigma_hat(2,2)/n) mu2_hat+sqrt((p*(n-1)/(n-p))*finv(1-alpha,p,n-p))*sqrt(Sigma_hat(2,2)/n)];
    if (mu1>=mu1_sim_CI(1) && mu1<=mu1_sim_CI(2) && mu2>=mu2_sim_CI(1) && mu2<=mu2_sim_CI(2))
        mu_in_sim_CI = mu_in_sim_CI + 1;
    end
    mu1_marg_CI = [mu1_hat-tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(1,1)/n) mu1_hat+tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(1,1)/n)];
    mu2_marg_CI = [mu2_hat-tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(2,2)/n) mu2_hat+tinv(1-alpha/2,n-1)*sqrt(Sigma_hat(2,2)/n)];
    if (mu1>=mu1_marg_CI(1) && mu1<=mu1_marg_CI(2) && mu2>=mu2_marg_CI(1) && mu2<=mu2_marg_CI(2))
        mu_in_marg_CI = mu_in_marg_CI + 1;
    end
    mu1_bonf_CI = [mu1_hat-tinv(1-alpha/(2*p),n-1)*sqrt(Sigma_hat(1,1)/n) mu1_hat+tinv(1-alpha/(2*p),n-1)*sqrt(Sigma_hat(1,1)/n)];
    mu2_bonf_CI = [mu2_hat-tinv(1-alpha/(2*p),n-1)*sqrt(Sigma_hat(2,2)/n) mu2_hat+tinv(1-alpha/(2*p),n-1)*sqrt(Sigma_hat(2,2)/n)];
    if (mu1>=mu1_bonf_CI(1) && mu1<=mu1_bonf_CI(2) && mu2>=mu2_bonf_CI(1) && mu2<=mu2_bonf_CI(2))
        mu_in_bonf_CI = mu_in_bonf_CI + 1;
    end
end
subplot(4,2,8)
axis off
text(-0.1,0.6,['Test of CR and CI''s with N = ' num2str(N) ' iterations each using n = ' num2str(n) ' observations:'],'Fontsize',14)
text(-0.1,0.4,['Desired confidence: 100(1-{\alpha})% = ' num2str(100*(1-alpha),3) ' %'],'Fontsize',14)
text(0.65,0.4,'Empirical values achieved:','Fontsize',14)
text(0.65,0.25,'{\mu} {\in} CR:','Fontsize',14)
text(0.89,0.25,[num2str(100*mu_in_CR/N,3) ' %'],'Fontsize',14)
text(1.04,0.25,'(black)','Fontsize',14,'Color','k')
text(0.65,0.1,'{\mu} {\in} simult.CI:','Fontsize',14)
text(0.89,0.1,[num2str(100*mu_in_sim_CI/N,3) ' %'],'Fontsize',14)
text(1.04,0.1,'(red)','Fontsize',14,'Color','r')
text(0.65,-0.05,'{\mu} {\in} marg.CI:','Fontsize',14)
text(0.89,-0.05,[num2str(100*mu_in_marg_CI/N,3) ' %'],'Fontsize',14)
text(1.04,-0.05,'(blue)','Fontsize',14,'Color','b')
text(0.65,-0.2,'{\mu} {\in} Bonf.CI:','Fontsize',14)
text(0.89,-0.2,[num2str(100*mu_in_bonf_CI/N,3) ' %'],'Fontsize',14)
text(1.04,-0.2,'(green)','Fontsize',14,'Color','g')
%--------------------------------------------------------------------------
% Text for figure
%--------------------------------------------------------------------------
subplot(4,2,1)
axis off
text(-0.3,1,'Confidence region (CR) and confidence intervals (CI) for {\mu} = ({\mu}_1,{\mu}_2)^T','Fontsize',20)
text(-0.3,0.75,'of an assumed bivariate normal, (X_1,X_2) ~ N_2({\mu},{\Sigma})   ({\mu},{\Sigma} unknown)','Fontsize',20)
text(-0.3,0.5,'Level of confidence:  100(1-{\alpha}) %','Fontsize',16)
text(-0.3,0.2,'CR (two-dimens.):     (x-{\mu})^T*S^-^1*(x-{\mu})  {\leq}  p(n-1)/n(n-p) * F_p_,_n_-_p({\alpha})','Fontsize',16)
text(-0.095,0.37,'                 _              _','Fontsize',16)
text(-0.3,-0.05,'Simultaneous CI''s:    ( x_i   {\pm}  sqrt [ p(n-1)/(n-p) * F_p_,_n_-_p({\alpha}) ] * sqrt [ s_i_i/n ] ),   i=1,2','Fontsize',16)
text(-0.093,0.13,'                 _','Fontsize',16)
text(-0.3,-0.3,'Marginal CI''s:            ( x_i   {\pm}  t_n_-_1({\alpha}/2) * sqrt [ s_i_i/n ] ),   i=1,2','Fontsize',16)
text(-0.093,-0.12,'                 _','Fontsize',16)
text(-0.3,-0.55,'Bonferroni CI''s:         ( x_i   {\pm}  t_n_-_1({\alpha}/2p) * sqrt [ s_i_i/n ] ),   i=1,2','Fontsize',16)
text(-0.093,-0.37,'                 _','Fontsize',16)
%------------------------------------------------------------------------------------------

