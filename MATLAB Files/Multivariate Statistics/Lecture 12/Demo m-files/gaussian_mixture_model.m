clc;clear;close all
format compact

figure
subplot(7,4,1:4)
axis off
text(0,1.25,'Finite Gaussian Mixture Model (GMM) for X:   f_X(x) = {\Sigma}_k_=_1_._._K[p_kf_Y_k(x)],     x = [x_1 ... x_p]^T','Fontsize',24)
text(0.43,0.55,'where    Y_k ~ N_p({\mu}_k,{\Sigma}_k),    {\Sigma}_k_=_1_._._K[p_k] = 1,    p_k {\geq} 0','Fontsize',24)
text(0,-0.1,'(note that f_X(x) is a valid p-dimensional pdf, f_X(x) {\geq} 0, {\int}f_X(x)dx = 1,  BUT X is NOT MVN (as apparent below)  and  X {\neq} {\Sigma}_k_=_1_._._K[p_kY_k])','Fontsize',14)
text(0,-0.65,'Examples for K = 2, p = 2','Fontsize',24)

%----------------------------
subplot(7,4,[9 13])
p1 = 0.5;
mu1 = [4 0];
Sigma1 = [4 0; 0 4];
%---------------------------
p2 = 1 - p1;
mu2 = [-4 0];
Sigma2 = [4 0; 0 4];
%----------------------------
x1 = -10:0.5:10;
x2 = -10:0.5:10;
[X1,X2] = meshgrid(x1,x2);
F = p1*mvnpdf([X1(:) X2(:)],mu1,Sigma1) + p2*mvnpdf([X1(:) X2(:)],mu2,Sigma2);
F = reshape(F,length(x2),length(x1));
%------------------------------
surf(x1,x2,F);
caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
xlabel('x_1','Fontsize',14); ylabel('x_2','Fontsize',14)
%--------------------------------
subplot(7,4,[10 14])
contour(x1,x2,F),axis square,axis([-10 10 -10 10]),xlabel('x_1','Fontsize',14); ylabel('x_2','Fontsize',14)
%-------------------------------  
subplot(7,4,5:6)
axis off
text(0,0.15,'pdf and contourplot for','Fontsize',14)
text(0.38,0.15,['p_1 = ' num2str(p1) ', {\mu}_1 = ' mat2str(mu1) ', {\Sigma}_1 = ' mat2str(Sigma1)],'Fontsize',14)
text(0.38,-0.15,['p_2 = ' num2str(p2) ', {\mu}_2 = ' mat2str(mu2) ', {\Sigma}_2 = ' mat2str(Sigma2)],'Fontsize',14)
%-------------------------------

%-------------------------------
subplot(7,4,[11 15])
p1 = 0.4;
mu1 = [4 0];
Sigma1 = [4 0; 0 4];
%---------------------------
p2 = 1 - p1;
mu2 = [-4 0];
Sigma2 = [4 0; 0 4];
%----------------------------
x1 = -10:0.5:10;
x2 = -10:0.5:10;
[X1,X2] = meshgrid(x1,x2);
F = p1*mvnpdf([X1(:) X2(:)],mu1,Sigma1) + p2*mvnpdf([X1(:) X2(:)],mu2,Sigma2);
F = reshape(F,length(x2),length(x1));
%------------------------------
surf(x1,x2,F);
caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
xlabel('x_1','Fontsize',14); ylabel('x_2','Fontsize',14)
%------------------------------
subplot(7,4,[12 16])
contour(x1,x2,F),axis square,axis([-10 10 -10 10]),xlabel('x_1','Fontsize',14); ylabel('x_2','Fontsize',14)
%------------------------------- 
subplot(7,4,7:8)
axis off
text(0,0.15,'pdf and contourplot for','Fontsize',14)
text(0.38,0.15,['p_1 = ' num2str(p1) ', {\mu}_1 = ' mat2str(mu1) ', {\Sigma}_1 = ' mat2str(Sigma1)],'Fontsize',14)
text(0.38,-0.15,['p_2 = ' num2str(p2) ', {\mu}_2 = ' mat2str(mu2) ', {\Sigma}_2 = ' mat2str(Sigma2)],'Fontsize',14)
%-------------------------------

%-------------------------------
subplot(7,4,[21 25])
p1 = 0.5;
mu1 = [2.5 0];
Sigma1 = [4 0; 0 1];
%---------------------------
p2 = 1 - p1;
mu2 = [-2.5 0];
Sigma2 = [1 0; 0 4];
%----------------------------
x1 = -10:0.5:10;
x2 = -10:0.5:10;
[X1,X2] = meshgrid(x1,x2);
F = p1*mvnpdf([X1(:) X2(:)],mu1,Sigma1) + p2*mvnpdf([X1(:) X2(:)],mu2,Sigma2);
F = reshape(F,length(x2),length(x1));
%------------------------------
surf(x1,x2,F);
caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
xlabel('x_1','Fontsize',14); ylabel('x_2','Fontsize',14)
%--------------------------------
subplot(7,4,[22 26])
contour(x1,x2,F),axis square,axis([-10 10 -10 10]),xlabel('x_1','Fontsize',14); ylabel('x_2','Fontsize',14)
%-------------------------------
subplot(7,4,17:18)
axis off
text(0,0.15,'pdf and contourplot for','Fontsize',14)
text(0.38,0.15,['p_1 = ' num2str(p1) ', {\mu}_1 = ' mat2str(mu1) ', {\Sigma}_1 = ' mat2str(Sigma1)],'Fontsize',14)
text(0.38,-0.15,['p_2 = ' num2str(p2) ', {\mu}_2 = ' mat2str(mu2) ', {\Sigma}_2 = ' mat2str(Sigma2)],'Fontsize',14)
%-------------------------------

%-------------------------------
subplot(7,4,[23 27])
p1 = 0.75;
mu1 = [0 -2];
Sigma1 = [4 -1; -1 4];
%---------------------------
p2 = 1 - p1;
mu2 = [-3 3];
Sigma2 = [1.4 1.3; 1.3 3];
%----------------------------
x1 = -10:0.5:10;
x2 = -10:0.5:10;
[X1,X2] = meshgrid(x1,x2);
F = p1*mvnpdf([X1(:) X2(:)],mu1,Sigma1) + p2*mvnpdf([X1(:) X2(:)],mu2,Sigma2);
F = reshape(F,length(x2),length(x1));
%------------------------------
surf(x1,x2,F);
caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
xlabel('x_1','Fontsize',14); ylabel('x_2','Fontsize',14)
%--------------------------------
subplot(7,4,[24 28])
contour(x1,x2,F),axis square,axis([-10 10 -10 10]),xlabel('x_1','Fontsize',14); ylabel('x_2','Fontsize',14)
%-------------------------------
subplot(7,4,19:20)
axis off
text(0,0.15,'pdf and contourplot for','Fontsize',14)
text(0.38,0.15,['p_1 = ' num2str(p1) ', {\mu}_1 = ' mat2str(mu1) ', {\Sigma}_1 = ' mat2str(Sigma1)],'Fontsize',14)
text(0.38,-0.15,['p_2 = ' num2str(p2) ', {\mu}_2 = ' mat2str(mu2) ', {\Sigma}_2 = ' mat2str(Sigma2)],'Fontsize',14)
%-------------------------------

