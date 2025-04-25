clc;clear;close all

mu1 = 0;
mu2 = 0;
sigma1 = 1;
sigma2 = 1;
var1 = sigma1^2;
var2 = sigma2^2;
rho = 0.7;

mu = [mu1 mu2];
Sigma = [var1 rho*sqrt(var1*var2); rho*sqrt(var1*var2) var2];

x1 = mu1-3*sigma1:6*sigma1/50:mu1+3*sigma1; 
x2 = mu2-3*sigma2:6*sigma2/50:mu2+3*sigma2;
[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],mu,Sigma);
F = reshape(F,length(x2),length(x1));

figure
subplot(2,2,3)
surf(x1,x2,F);
caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
axis([-3 3 -3 3 0 .4])
xlabel('x_1','Fontsize',14); ylabel('x_2','Fontsize',14);title('3D plot of probability density function   f(x_1, x_2)','Fontsize',14)

subplot(2,2,2)
maxvar = max(var1,var2);
contour(x1,x2,F),axis square,title({'Probability density contour plot for f(x_1,x_2)';'and ellipse eigenvector axes  {\surd}{\lambda}_i * e_i  (shown for c^2=1) '},'Fontsize',14),grid,...
    axis([mu1-3*sqrt(maxvar) mu1+3*sqrt(maxvar) mu2-3*sqrt(maxvar) mu2+3*sqrt(maxvar)]),xlabel('x_1','Fontsize',14); ylabel('x_2','Fontsize',14)
hold on
[V, lambda] = eig(Sigma);
line([mu1  mu1+sqrt(lambda(1,1))*V(1,1)],[mu2 mu2+sqrt(lambda(1,1))*V(2,1)],'LineWidth',4,'Color','k')
line([mu1  mu1+sqrt(lambda(2,2))*V(1,2)],[mu2 mu2+sqrt(lambda(2,2))*V(2,2)],'LineWidth',4,'Color','k')
hold off

subplot(224)
N = 1000;
obs = mvnrnd(mu,Sigma,N);
alpha = [0.5 0.25 0.05];
plot(obs(:,1),obs(:,2),'+'),grid,axis([mu1-3*sqrt(maxvar) mu1+3*sqrt(maxvar) mu2-3*sqrt(maxvar) mu2+3*sqrt(maxvar)]),...
    xlabel('x_1','Fontsize',14); ylabel('x_2','Fontsize',14),...
    title({['N = ' num2str(N) ' observations scatter plot for (X_1, X_2) and'];['(1-{\alpha})*100% prediction ellipses for  {\alpha}  =  [' num2str(fliplr(alpha)) ']']},'Fontsize',14)
hold on

df = 2;
c = chi2inv(1-alpha,df);
C = (1/(2*pi*sqrt(det(Sigma)))).*exp(-c/2);
contour(x1,x2,F,C,'Color','k','LineWidth',2),axis square
hold off

subplot(221)
axis off
text(0.1,1.2,'Bivariate normal distribution parameters and "geometry"','Fontsize',14)
text(0.1,1.09,'(X_1,X_2) ~ N({\mu},{\Sigma})','Fontsize',12)
text(0.1,1,['{\mu} = [' num2str(mu1) ', ' num2str(mu2) ']'],'Fontsize',12)
text(0.3,1,['[{\sigma}_1, {\sigma}_2 ] = [' num2str(sigma1) ', ' num2str(sigma2) ']'],'Fontsize',12) 
text(0.6,1,['{\rho} = ' num2str(rho)],'Fontsize',12)
text(0.1,0.92,['{\Sigma} = [' num2str(Sigma(1,1)) '  ' num2str(Sigma(1,2)) '; ' num2str(Sigma(2,1)) '  ' num2str(Sigma(2,2))  ']'],'Fontsize',12) 
text(0.1,0.83,['total variation                =   {\Sigma}_i {\lambda}_i   =   {\sigma}_1^2 + {\sigma}_2^2   =   ' num2str(Sigma(1,1)+Sigma(2,2))],'Fontsize',12)
text(0.1,0.73,['generalized variance   =   | {\Sigma} |    =   {\Pi}_i {\lambda}_i  =   {\sigma}_1^2 * {\sigma}_2^2 * (1 - {\rho}^2)   =   ' num2str(det(Sigma))],'Fontsize',12)
text(0.1,0.65,'------------------------------------------------------------------------------------------------------------------------------------')
text(0.1,0.58,'constant density contour ellipses   (x-{\mu})^T*{\Sigma}^-^1*(x-{\mu}) = c^2','Fontsize',12)
text(0.1,0.46,['eigenvalues and vectors for {\Sigma}        {\lambda}_1 = ' num2str(lambda(1,1)) '       e_1 = [' num2str(V(1,1),3) '  ' num2str(V(2,1),3) ']'],'Fontsize',12)  
text(0.1,0.38,['                                                           {\lambda}_2 = ' num2str(lambda(2,2)) '       e_2 = [' num2str(V(1,2),3) '  ' num2str(V(2,2),3) ']'],'Fontsize',12)
text(0.1,0.32,'------------------------------------------------------------------------------------------------------------------------------------')
text(0.1,0.23,'Mahalanobis (statistical) distance   (x-{\mu})^T*{\Sigma}^-^1*(x-{\mu}) ~ {\chi}^2_2','Fontsize',12)
text(0.1,0.16,'------------------------------------------------------------------------------------------------------------------------------------')
text(0.1,0.08,'therefore  c^2 = {\chi}^2_2_,_{\alpha}  gives (1-{\alpha})*100% prediction ellipses for observations,  x_i','Fontsize',12)
text(0.1,-0.02,['ellipse half-length axes    {\surd}{\lambda}_i * {\surd}{\chi}^2_2_,_{\alpha} * e_i   =    [' num2str(sqrt(lambda(1,1)),3) ' * e_1,   '  num2str(sqrt(lambda(2,2)),3) ' *e_2] * {\surd}{\chi}^2_2_,_{\alpha}'],'Fontsize',12)
text(0.1,-0.12,['ellipse areas  A_{\alpha} =  {\pi} * {\chi}^2_2_,_{\alpha} * {\surd}|{\Sigma}|   =   [' num2str(pi*fliplr(c)*sqrt(det(Sigma)),3) ']'],'Fontsize',12)
text(0.1,-0.19,['                                              for   {\alpha}    =   [' num2str(fliplr(alpha)) ']'],'Fontsize',12)
text(0.1,-0.24,'------------------------------------------------------------------------------------------------------------------------------------')

