clc
clear
close all

syms k
sys = tf([1],[1 1 0]);
K = tf([1 2],[1 10]);

%figure
%rlocus(sys*K)

%1+k*K*sys = 1+k*G = 0

k = 1.6

sysP = sys*k*K/(1+sys*k*K);

%step(sysP)
%pole(sysP)
%zero(sysP)

syms s m g l I b theta

C = [1 0];

M = [s -1;
     (-m*g*l*cos(theta))/I s+(b/I)]
 
M = simplify(inv(M))

B = [0 0;
     0 1/I];
D = 0;

H = C * M * B + D

H = simplify(H)

pretty(H(2))
