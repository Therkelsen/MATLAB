clc; clear; close all;
format compact;

syms m_dot u q1 q2 q3 omega omega_dot m J
eq1 = m_dot == u - q1 * omega * m
eq2 = J * omega_dot == q2 * m - q3 * omega

% Define parameters and equilibrium values
params = [q1, q2, q3, J];  % parameters in equations
values = [1, 10, 0.1, 40];  % actual values for q1, q2, q3, J

%1. Characterize (linear/nonlinear, order, control input vector dimension)
% and give a state-space representation of system

% Linearity:
% Note: Look for products, powers, trigonometric functions,
% or other nonlinear operations of state variables or inputs
% Answer: Both equations have products, showing nonlinearity.

% Order:
% Note: Count the number of first-order differential equations (m_dot and
% omega_dot are given here).
% Each state variable (m and ω) contributes to the order.
% Answer: Two coupled first-order equations,
% meaning this is a second order system.

% Control Input Vector Dimension:
% Note: Look at the terms directly involving u. Check how many independent
% inputs influence the system.
% Answer: For the system: The control input u appears only in the first
% equation (m_dot = u - q_1*omega*m). There is only one input, with dim 1.

% State Space Representation:
% Note: State variables should be chosen to fully describe the dynamics
% of the system at any time. The state-space representation allows you to
% model the system's evolution using state variables.
% 
% Guidelines for selecting state variables:
%  - 1st-order ODE: Every differential equation that describes the time
%    evolution of a variable should be associated with a state variable.
%  - Independent Variables: Choose variables that describe the system’s
%    behavior. Typically, you choose the position, velocity, or any
%    relevant quantities that encapsulate the system's state.
%
% In your system, you have:
%  - m (mass flow rate) as a state variable: Since it appears in a
%    differential equation, it represents a dynamic quantity that
%    needs to be tracked over time.
%  - ω (angular velocity) as a state variable: Similarly, since it
%    is part of the second differential equation and directly affects
%    system behavior, it should also be a state variable.

m_dot = rhs(eq1)
omega_dot = rhs(eq2)/J

% Substitute the parameters here
m_dot = subs(m_dot, params, values)
omega_dot = subs(omega_dot, params, values)

X = [m; omega]
X_dot = [m_dot; omega_dot]

% 2. Simulink

% 3. Give the set of equilibrium points of system (1)-(2)
% and their associated control inputs.

% Infinite EQ points.

syms u q1 x1 x2
equ1 = u == q1*x1*x2;
equ2 = x2 == q2/q3*x1;
equ3 = x1 == q3/q2*x2;
omega = rhs(equ2);
m = rhs(equ3);
u = subs(rhs(equ1),x1, m);

omega = subs(omega, params, values)
m = subs(m, params, values)
u = subs(u, params, values)

% For omega_star = 200
omega_star_1 = 200
m_star_1 = double(subs(m, x2, omega_star_1))
u_star_1 = double(subs(u, x2, omega_star_1))

omega_star_2 = 300
m_star_2 = double(subs(m, x2, omega_star_2))
u_star_2 = double(subs(u, x2, omega_star_2))

syms m omega u
X_dot = subs(X_dot, {m, omega}, {x1, x2})

df1dx1 = diff(X_dot(1), x1);
df1dx2 = diff(X_dot(1), x2);
df2dx1 = diff(X_dot(2), x1);
df2dx2 = diff(X_dot(2), x2);

A = [df1dx1, df1dx2;
     df2dx1, df2dx2]

A = double(subs(A, {x1, x2}, {m_star_2, omega_star_2}))

B = diff(X_dot, u)

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

% Stable, but only barely.