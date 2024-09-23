clc; clear; close all;
format compact

% Exercise 1
s = tf('s')


theta = 160*(s + 2.5)*(s + 0.7);
delta = (s*s + 5*s + 40)*(s*s + 0.03*s + 0.06);

P = theta/delta

isstable(P)

figure(1)
bode(P)
grid on
% 
% figure(2)
% pzplot(P)
% grid on

% Exercise 2, 3, 4, 5
% Plus stuff in Simulink

pg = getPeakGain(P)

dg = dcgain(P)

Kp = 0.5;
Kd = 0.25;
Ki = 0.05;

% Exercise 6
C = Kp + Ki/s + Kd*s

H = C*P

H2 = feedback(C * P, 1)

% Exercise 7

figure(3)
bode(H)
grid on

% Exercise 8
