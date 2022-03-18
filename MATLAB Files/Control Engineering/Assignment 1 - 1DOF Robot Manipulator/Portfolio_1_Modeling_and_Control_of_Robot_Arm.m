
%% PID Controller Root Locus
clear
close all
G = tf(1,[1/3,0.1,-2.455]);
s = tf('s');

sigma = -0.28;
zeta = 0.3495;
omega_n = 0.225;

N = 10;
Kp = 1;
Ti = 1;
Td = 1/5;

%K = Kp*(1 + 1/(Ti*s) + Td*s);
K = Kp*(1 + 1/(s*Ti) + (s*Td)/(1+s*(Td/N)));

L = K*G;

zero(L/(1+L))

figure(1)
rlocus(L)

figure(2)
step(K)
%%
% plot requirenments %
hold on
plot([sigma sigma],[-5 5], '--', 'color', 'r')                      % sigma
plot([0 5*(-zeta)],[0 5*cos(asin(-zeta))], '--', 'color', 'r')      % zeta
plot([0 5*(-zeta)],[0 -5*cos(asin(-zeta))], '--', 'color', 'r')     % zeta
th = linspace( pi/2, -pi/2, 100);                                   % omega_n
R = omega_n;
x = -R*cos(th) + 0;
y = -R*sin(th) + 0;
plot(x,y,'--', 'color', 'r');
axis equal;


