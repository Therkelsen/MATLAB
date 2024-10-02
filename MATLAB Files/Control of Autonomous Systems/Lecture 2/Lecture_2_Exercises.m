clc; close all;
format compact

% Exercise 1.1, 1.2, 1.3

x_trace = squeeze(out.x); % transform 3D array into 2D array

close all % close all windows

% figure, plot3(x_trace(1,:),x_trace(2,:),x_trace(3,:))
% grid on
% xlabel('x')
% ylabel('y')
% zlabel('z')
% title('Evolution of state for Lorenz system')

% Exercise 2.1, 2.2, 2.3

% Solution to Exercise 2 of Lecture 2: ODE to SS rep
% Jerome Jouffroy, September 2024
% the ODE is y_dddot - 4 y_ddot + 7 y_dot - 2 y = 3 u
% LHS is denominator, RHS is numerator

s = tf('s') % definition of s as Laplace variable

H = 3/(s^3-4*s^2+7*s-2) % definition of transfer function

[num_H,den_H] = tfdata(H,'v') % extraction of numerator and denominator of TF function H

[A,B,C,D] = tf2ss(num_H,den_H) % transform into state-space representation

% Exercise 3.1
m = 50;
l = 1;
g = 9.8;
d = 0.1;

% parameter vector
a = [m;
     g;
     l;
     d]

% equilibrium points

xstar = [pi/4;
         0]

ustar = m*g*l*sin(xstar(1))

% initial condition
x0 = xstar;

% linear approximation of pendulum around xstar (and associated ustar)
A = [ 0 1 ; -g/l*cos(pi/4) -d/(m*l^2) ]
B = [ 0 ; 1/(m*l^2) ]