clc; clear; close all;
format compact;

syms Fc gamma I I_dot J Ke Km Kd L omega omega_dot Rm u
eq1 = u == L*I_dot + Rm*I + Ke*omega
eq2 = Km * I == J*omega_dot + Kd*omega + Fc*gamma

X = [I; omega]

I_dot = simplify(solve(rhs(eq1 - u), I_dot))
I_dot = expand(I_dot)

omega_dot = simplify(solve(eq2, omega_dot))
omega_dot = expand(omega_dot)

X_dot = [I_dot; omega_dot]

params = [Fc];
values = [0];

syms x1 x2

states = [I, omega];
ss_states = [x1, x2;]

X = subs(X, states, ss_states)
X_dot = subs(X_dot, params, values)
X_dot = subs(X_dot, states, ss_states)

u_eq = solve(X_dot(1) == 0, u);
disp('u_eq:')
pretty(u_eq)

I_eq = solve(X_dot(2) == 0, x1);
disp('I_eq:')
pretty(I_eq)

omega_eq = solve(X_dot(2) == 0, x2);
disp('omega_eq:')
pretty(omega_eq)

u_eq = subs(u_eq, x1, I_eq);
disp('x1 subbed u_eq:')
pretty(u_eq)

params = [Fc, J, L, Kd, Ke, Km, Rm];
values = [0, 0.01, 1, 0.001, 1, 1, 10];

X_dot = subs(X_dot, params, values)

% Stability

df1dx1 = diff(X_dot(1), x1)
df1dx2 = diff(X_dot(1), x2)
df2dx1 = diff(X_dot(2), x1)
df2dx2 = diff(X_dot(2), x2)

A = [df1dx1, df1dx2;
     df2dx1, df2dx2]
A = double(A);
% A = double(subs(A, {x1, x2}, {I_star_1, omega_star_1}))

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

% Controllability matrix
n = length(X);
C = B;  % Initialize with B
for i = 1:n-1
    C = [C, A^i * B];
end

disp('B:')
disp(B)

disp('A*B')
disp(A*B)

% Check rank of the controllability matrix
rankC = rank(C);  % Symbolic rank

% Output
disp('Controllability Matrix:');
disp(C);
disp(['Rank of Controllability Matrix: ', num2str(rankC)]);
if rankC == n
    disp('The system is controllable.')
else
    disp('The system is not controllable.')
end

% Q5 - Simulink Implementation
u_eq = subs(u_eq, params, values);
I_eq = subs(I_eq, params, values);
omega_eq = subs(omega_eq, params, values);

omega_star_1 = 10
I_star_1 = subs(I_eq, x2, omega_star_1)
I_star_1 = double(I_star_1)
x_star_1 = [I_star_1, omega_star_1];
u_star_1 = subs(u_eq, ss_states, x_star_1)
u_star_1 = double(u_star_1)

% System stays at x_star_1, it is indeed a valid eq point.

% Q6 State Feedback Control
disp(' ')
disp(' ')
disp('State Feedback Control')
% CANNOT HAVE TWO POLES THAT ARE THE SAME
% Multiplicity of poles must not exceed the rank of the B matrix
% Perhaps its better to pick complex-conjugated poles then?
lambda_cl = [-5 + 1i, -5 - 1i]

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

x_star_2 = [0, 10]

randseed = randi([])