clc
clear
close all

%% Model Parameters
m = 1; % Mass of link [kg]
l = 0.5; % Distance to CoM of link [m]
b = 0.1; % Joint friction coefficient [Nm/(rad/s)]
g = 9.82; % Gravitational acceleration [m/s^2]

%% Root Locus Plot
clear
%H = (s+1)/(s^2+2*s+2);
G = tf(1,[1/3,0.1,-2.455]);
s = tf('s');

sigma = -0.28;
zeta = 0.3495;
omega_n = 0.225;

N = 10;
Kp = 0.97;
Ti = 1;
Td = 1/5;

%K = Kp*(1 + 1/(Ti*s) + Td*s);
K = Kp*(1 + 1/(Ti*s) + (Td*s)/(1+s*Td/N));

zero((K*G)/(1+K*G))

H = (K*G)/(1+K*G);

[R,K] = rlocus(H);
polesH = pole(H);
zerosH = zero(H);

rlocus((1 + 1/(Ti*s) + (Td*s)/(1+s*Td/N)))

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
%% Step response
step(H, 30)

%% System Response
t = linspace(0,10);
u = sin(t);
y = lsim(H ,u ,t);
%y = 0.5*sin(t+1);

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