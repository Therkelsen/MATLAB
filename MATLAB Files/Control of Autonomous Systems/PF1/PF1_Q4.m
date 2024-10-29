clc; clear; close all;
format compact

syms P1 P2 u Ts

X = [P1;
     P2]

Xdot = [X(1)*X(2)-X(1);
        -4/3*X(1)*X(2)+2/3*X(2)+u]

d_X = [X(1) + Ts*(Xdot(1));
       X(2) + Ts*(Xdot(2))]

syms x1 x2

Xdot = subs(subs(Xdot, P1, x1), P2, x2)

dfdx = [diff(Xdot(1),x1) diff(Xdot(1),x2);
        diff(Xdot(2),x1), diff(Xdot(2), x2)]

matlabFunction(dfdx, 'File', 'computeA', 'Vars', {x1, x2});

% Q7
xeq1 = [0;
        0];

A1 = double(subs(subs(dfdx, x1, xeq1(1)), x2, xeq1(2)))

xeq2 = [1/2;
        1];

A2 = double(subs(subs(dfdx, x1, xeq2(1)), x2, xeq2(2)))

B = double([diff(Xdot(1), u);
            diff(Xdot(2), u)])

e = eig(A2)
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

% Q8
% Desired poles for the closed-loop system
disp(' ')
disp(' ')
disp('State Feedback Control')
lambda_cl = [-1, -2]

K = place(A2, B, lambda_cl)

A_cl = A2 - B * K

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

% State feedback
xeq3 = [5;
        1];

A3 = double(subs(subs(dfdx, x1, xeq3(1)), x2, xeq3(2)))

B = double([diff(Xdot(1), u);
            diff(Xdot(2), u)])

e = eig(A3)
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

% Controller
% Desired poles for the closed-loop system
disp(' ')
disp(' ')
disp('State Feedback Control')
lambda_cl = [-1, -2]

K2 = place(A3, B, lambda_cl)

A_cl = A3 - B * K2

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

% Integral Feedback Control
C = [1 0];

Ki = -10;

