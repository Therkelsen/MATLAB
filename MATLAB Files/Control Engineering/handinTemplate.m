clc
clear
close all

%% Model Parameters
m = 1; % Mass of link [kg]
l = 0.5; % Distance to CoM of link [m]
b = 0.1; % Joint friction coefficient [Nm/(rad/s)]
g = 9.82; % Gravitational acceleration [m/s^2]

%% Root Locus Plot
s = tf('s');
H = (s+1)/(s^2+2*s+2);
[R,K] = rlocus(H);
polesH = pole(H);
zerosH = zero(H);


font_size=10;
width = 10;         
height = 5;
SCREEN_LEFT = 15;
SCREEN_RIGHT = 10;
ADD = 0;
figure
hold on
set(gcf,'Units','centimeters')
set(gcf,'Position',[SCREEN_LEFT  SCREEN_RIGHT  width+4*ADD  height+2*ADD])
set(gcf,'Color','w')
hold on
ax = gca;

plot([-1e6 1e6],[0 0],':k')
plot([0 0],[-1e6 1e6],':k')
ax.ColorOrderIndex = 1;
plot(R')
ax.ColorOrderIndex = 1;
plot(complex(real(polesH),imag(polesH)),'x')
ax.ColorOrderIndex = 1;
plot(complex(real(zerosH),imag(zerosH)),'o')
set(gca,'FontName','Times New Roman','FontSize',font_size);
xlabel('Re$(s)$','Interpreter','latex','FontSize',font_size)
ylabel('Im$(s)$','Interpreter','latex','FontSize',font_size)
xlim([-4 4])
ylim([-2 2])

box on

%% System Response
t = linspace(0,10);
u = sin(t);
y = 0.5*sin(t+1);

font_size=10;
width = 10;         
height = 5;
SCREEN_LEFT = 15;
SCREEN_RIGHT = 10;
ADD = 0;
figure
hold on
set(gcf,'Units','centimeters')
set(gcf,'Position',[SCREEN_LEFT  SCREEN_RIGHT  width+4*ADD  height+2*ADD])
set(gcf,'Color','w')
subplot(2,1,1)
hold on
plot(t,y)
set(gca,'FontName','Times New Roman','FontSize',font_size);
set(gca,'XTickLabel','');
ylabel('Output $y$','Interpreter','latex','FontSize',font_size)
box on
subplot(2,1,2)
hold on
plot(t,u)
set(gca,'FontName','Times New Roman','FontSize',font_size);
ylabel('Input $u$','Interpreter','latex','FontSize',font_size)
xlabel('Time $t$ [s]','Interpreter','latex','FontSize',font_size)


box on
