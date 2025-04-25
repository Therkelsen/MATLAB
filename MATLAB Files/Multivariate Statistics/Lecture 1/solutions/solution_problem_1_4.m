clc;clear;close all
format compact

n = 100;
load('M:\Statistik\MULTISTAT\Lektion 01  Multivariate Data\Course material\dataset_problem_1_4.mat');

figure
subplot(331)
axis off
text(-0.5,1,'Kernel Density Estimation of unknown pdf','Fontsize',20)
text(-0.5,0.75,'Random sample:             x_1, x_2, ..., x_n','Fontsize',16)
text(-0.5,0.6,'Kernel:                             K(x),  where {\int}K(x)dx = 1','Fontsize',16)
text(-0.5,0.4,'Kernel width parameter:   h','Fontsize',16)
text(-0.5,0.1,'pdf estimate:   f_K_D_E(x) = (1/nh){\cdot}{\Sigma}_i_=_1_._._nK[(x-x_i)/h]','Fontsize',20)

subplot(332)
plot(xplot,f,'r','Linewidth',3),grid,title('true pdf, f(x)  (gaussian mixture model)','Fontsize',14),ylim([0 0.4])

subplot(333)
hist(x,30)
title(['histogram of n = ' num2str(n) ' random sample'],'Fontsize',14)
 
kernel1 = 'normal';

h11 = 0.14;
fKDE = ksdensity(x,xplot,'kernel',kernel1,'width',h11,'function','pdf');
AvRMSerror = sqrt(mean((f-fKDE).^2));
subplot(334)
plot(xplot,f,'r','Linewidth',2),grid,title(['f_K_D_E(x)  with K = ' kernel1 ', h = ' num2str(h11)],'Fontsize',14)
hold on
plot(xplot,fKDE,'k','Linewidth',2),ylim([0 0.4])
text(2,0.35,['Average RMS error = ' num2str(AvRMSerror)],'Fontsize',12)
hold off

h12 = 0.26;
fKDE = ksdensity(x,xplot,'kernel',kernel1,'width',h12,'function','pdf');
AvRMSerror = sqrt(mean((f-fKDE).^2));
subplot(335)
plot(xplot,f,'r','Linewidth',2),grid,title(['f_K_D_E(x)  with K = ' kernel1 ', h = ' num2str(h12)],'Fontsize',14)
hold on
plot(xplot,fKDE,'k','Linewidth',2),ylim([0 0.4])
text(2,0.35,['Average RMS error = ' num2str(AvRMSerror)],'Fontsize',12)
hold off

h13 = 0.38;
fKDE = ksdensity(x,xplot,'kernel',kernel1,'width',h13,'function','pdf');
AvRMSerror = sqrt(mean((f-fKDE).^2));
subplot(336)
plot(xplot,f,'r','Linewidth',2),grid,title(['f_K_D_E(x)  with K = ' kernel1 ', h = ' num2str(h13)],'Fontsize',14)
hold on
plot(xplot,fKDE,'k','Linewidth',2),ylim([0 0.4])
text(2,0.35,['Average RMS error = ' num2str(AvRMSerror)],'Fontsize',12)
hold off

kernel2 = 'box';

h21 = 0.14;
fKDE = ksdensity(x,xplot,'kernel',kernel2,'width',h21,'function','pdf');
AvRMSerror = sqrt(mean((f-fKDE).^2));
subplot(337)
plot(xplot,f,'r','Linewidth',2),grid,title(['f_K_D_E(x)  with K = ' kernel2 ', h = ' num2str(h21)],'Fontsize',14)
hold on
plot(xplot,fKDE,'b','Linewidth',2),ylim([0 0.4])
text(2,0.35,['Average RMS error = ' num2str(AvRMSerror)],'Fontsize',12,'Color','b')
hold off

h22 = 0.26;
fKDE = ksdensity(x,xplot,'kernel',kernel2,'width',h22,'function','pdf');
AvRMSerror = sqrt(mean((f-fKDE).^2));
subplot(338)
plot(xplot,f,'r','Linewidth',2),grid,title(['f_K_D_E(x)  with K = ' kernel2 ', h = ' num2str(h22)],'Fontsize',14)
hold on
plot(xplot,fKDE,'b','Linewidth',2),ylim([0 0.4])
text(2,0.35,['Average RMS error = ' num2str(AvRMSerror)],'Fontsize',12,'Color','b')
hold off

h23 = 0.38;
fKDE = ksdensity(x,xplot,'kernel',kernel2,'width',h23,'function','pdf');
AvRMSerror = sqrt(mean((f-fKDE).^2));
subplot(339)
plot(xplot,f,'r','Linewidth',2),grid,title(['f_K_D_E(x)  with K = ' kernel2 ', h = ' num2str(h23)],'Fontsize',14)
hold on
plot(xplot,fKDE,'b','Linewidth',2),ylim([0 0.4])
text(2,0.35,['Average RMS error = ' num2str(AvRMSerror)],'Fontsize',12,'Color','b')
hold off

%---------------------------------------------------------------------------
