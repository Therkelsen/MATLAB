%% Simulink konstanter
linear_num = [1];
linear_denom = [1/12 0.1 -9.82/4];

Ti = 1;
Td = 1/7;
Kp = 30;

g = 9.82;
l = 0.5;
m = 1;
b = 0.1;
I = 1/12;

y_bar = pi/3;
u_bar = -(l*g*m*sqrt(3)) / 2;