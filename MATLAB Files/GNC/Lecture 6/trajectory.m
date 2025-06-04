clc; clear; close all;
format compact;

syms x y theta F1 F2
% Roomba states
X = [x; y; theta]

% Flatness variables chosen
F = [X(1);
     X(2)]

% States as a function of F and its' derivatives
X = subs(X, [x, y, theta], [F1, F2, atan(diff(F(1))/diff(F(2)))])


T = 10 % seconds

x0 = [0, 0, pi/4] % m, m, rad
xT = [1, 5, 3*pi/4] % m, m, rad

u10 = 1;
u1T = 1;

v0 = [u10*cos(x0(3)), u10*sin(x0(3)), 0] % m/s, m/s, rad/s
vT = [u10*cos(xT(3)), u10*sin(xT(3)), 0] % m/s, m/s, rad/s

Xpoly = [x0(1) x0(2) x0(3);
         v0(1) v0(2) v0(3);
         xT(1) xT(2) xT(3);
         vT(1) vT(2) vT(1)]

Bpoly = [b(0)';
         bp(0)';
         b(T)';
         bp(T)']

alpha = inv(Bpoly)*Xpoly

function b = b(t)
    b = [1, t, t^2, t^3]';
end

function bp = bp(t)
    bp = [0, 1, 2*t, 3*t^2]';
end