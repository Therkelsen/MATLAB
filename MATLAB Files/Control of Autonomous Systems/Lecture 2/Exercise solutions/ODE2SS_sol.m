% Solution to Exercise 2 of Lecture 2: ODE to SS rep
% Jerome Jouffroy, September 2024
% the ODE is y_dddot - 4 y_ddot + 7 y_dot - 2 y = 3 u

s = tf('s') % definition of s as Laplace variable

H = 3/(s^3-4*s^2+7*s-2) % definition of transfer function

[num_H,den_H] = tfdata(H,'v') % extraction of numerator and denominator of TF function H

[A,B,C,D] = tf2ss(num_H,den_H) % transform into state-space representation