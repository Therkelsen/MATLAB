% Parameters and system definition for the Piper example
% Jerome Jouffroy, September 2024

num1 = [1 2.5];
den1 = [1 5 40];
sys1 = tf(num1,den1)

num2 = [1 0.7];
den2 = [1 0.03 0.06];
sys2 = tf(num2,den2)

% cascade of 2 systems
P = 160*sys1*sys2

% plot its Bode diagrams
figure, bode(P)
grid on

% calculate DC gain
P_DCgain = dcgain(P)

% extraction of numerator and denominator
[ numP , denP ] = tfdata(P,'v');

% PID tuning
Kp = 0.5;
Ki = 0.005;
Kd = 1;

% obtain the transfer function of the PID controller
s = tf('s');
C = Kp + Ki/s + Kd*s

% Transfer function of closed-loop system
H = feedback(C*P,1)
isstable(H) % check whether the CL system is stable
H_DCgain = dcgain(H) % compute the DC gain
% Plot the Bode diagrams
figure, bode(H)
grid on
title('Transfer function of Piper stabilized with PID controller')


