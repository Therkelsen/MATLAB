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

params = [Fc, J, L, Kd, Ke, Km, Rm];
values = [0, 0.01, 1, 0.001, 1, 1, 10];

syms x1 x2

states = [I, omega];
ss_states = [x1, x2;]

X_dot = subs(X_dot, params, values)
X_dot = subs(X_dot, states, ss_states)

u_eq = solve(X_dot(1) == 0, u)
I_eq = solve(X_dot(2) == 0, x1)
omega_eq = solve(X_dot(2) == 0, x2)
u_eq = subs(u_eq, x1, I_eq)

omega_star_1 = 10
I_star_1 = double(subs(I_eq, x2, omega_star_1))
x_star_1 = [omega_star_1, I_star_1];
u_star_1 = double(subs(u_eq, ss_states, x_star_1))

% % Compute A and B matrices symbolically
% A = jacobian(X_dot, X);  % Partial derivatives w.r.t states
% B = jacobian(X_dot, u);  % Partial derivatives w.r.t input
% 
% % Controllability matrix
% n = length(X);
% C = B;  % Initialize with B
% for i = 1:n-1
%     C = [C, A^i * B];
% end
% 
% % Check rank of the controllability matrix
% rankC = rank(C);  % Symbolic rank
% 
% % Output
% disp('Controllability Matrix:');
% disp(C);
% disp(['Rank of Controllability Matrix: ', num2str(rankC)]);
% if rankC == n
%     disp('The system is controllable.')
% else
%     disp('The system is not controllable.')
% end