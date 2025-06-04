clc; clear;

syms F1 F2 F3 F4 F5 F6 L tau1 tau2 tau3 tau4 tau5 tau6 Ix Iy Iz K d b m real

U1 = F1 + F2 + F3 + F4 + F5 + F6 % Altitude
U2 = (F2- F1 + (F3 + F6 - F4 - F5)/2)*L % Roll
U3 = (F3 + F5 - F4 - F6)*L * sqrt(3) / 2 % Pitch
U4 = tau2 + tau5 + tau6 - tau1 - tau3 - tau4 % Yaw

U = [U2 U3 U4 U1]'

B = [0 0 0 0;
     0 0 0 0;
     0 0 0 0;
     0 0 0 0;
     0 -L/Iy 0 L/Ix; % Pitch and Roll -> Yaw Rate. Element 2 should by Iy, since it's pitch?
     -L/Iy 0 L/Iy 0;
     d/(b*Iz) -d/(b*Iz) d/(b*Iz) -d/(b*Iz); % Z Rate. Could make sense that it has to compensate for g, so it's nonzero.
     1/m 1/m 1/m 1/m] % Z Accel. Could make sense that it has to compensate for g, so it's nonzero.

BU = simplify(B*U)

BU = subs(BU, [F1, F2, F3, F4, F5, F6], [1, 1, 1, 1, 1, 1]);
disp(expand(BU))