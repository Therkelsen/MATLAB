clc; clear; close all;
format compact;

syms d u lambda x1 x2

X = [x1;
     x2]

X_dot = [-d*x1*abs(x1)+x2;
         -lambda*x2+u];

params = [d, lambda];
values = [0.005, 10];

X_dot = subs(X_dot, params, values)

u_eq = solve(X_dot(2) == 0, u);
disp('u_eq:')
pretty(u_eq)

x2_eq = solve(X_dot(1) == 0, x2);
disp('x2_eq:')
pretty(x2_eq)


x1_star_1 = 20
x2_star_1 = double(subs(x2_eq, x1, x1_star_1))
u_star_1 = double(subs(u_eq, x2, x2_star_1))

x_star_1 = [x1_star_1, x2_star_1]