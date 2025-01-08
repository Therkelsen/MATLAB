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

equ_1 = [u_star_1, x_star_1]

disp(' ')
disp(' ')
disp('Linearization around EQ Point')

x1_star_2 = 30
x2_star_2 = double(subs(x2_eq, x1, x1_star_2))
u_star_2 = double(subs(u_eq, x2, x2_star_2))

x_star_2 = [x1_star_2, x2_star_2]

equ_2 = [u_star_2, x_star_2]

df1dx1 = -abs(x1)/100;
df1dx2 = diff(X_dot(1), x2)
df2dx1 = diff(X_dot(2), x1)
df2dx2 = diff(X_dot(2), x2)

A = [df1dx1, df1dx2;
     df2dx1, df2dx2]

A = double(subs(A, {x1, x2}, x_star_2))

B = double(diff(X_dot, u))

e = eig(A)
poles = real(e);
% Check the real parts of the eigenvalues
disp(" ")
disp('Real parts of the eigenvalues: ')
disp(poles)
disp(" ")
if all(poles < 0)
    disp('The system is stable.');
elseif any(poles > 0)
    disp('The system is unstable.');
else
    disp('The system is marginally stable.');
end

% Stable.

disp(' ')
disp(' ')
disp('State Feedback Control')
lambda_cl = [-1, -10]

K = place(A, B, lambda_cl)

A_cl = A - B * K

e_cl = eig(A_cl)
poles_cl = real(e_cl);
% Check the real parts of the eigenvalues
disp(" ")
disp('Real parts of the eigenvalues: ')
disp(poles_cl)
disp(" ")
if all(poles_cl < 0)
    disp('The system is stable.');
elseif any(poles_cl > 0)
    disp('The system is unstable.');
else
    disp('The system is marginally stable.');
end