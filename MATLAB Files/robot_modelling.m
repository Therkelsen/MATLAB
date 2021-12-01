% Rotation, Translation, Transformation & Center of mass for links
clc
clear

syms MMP01 MMP02 A01 A12 R01 R12 P01 P12 q1 q2 L1 L2 Lc1 Lc2;

% Rotation matrices
R01 = [cos(q1) -sin(q1) 0; sin(q1) cos(q1) 0; 0 0 1];
R12 = [cos(q2) -sin(q2) 0; sin(q2) cos(q2) 0; 0 0 1];

% Translation vectors
P01 = [L1*cos(q1) L1*sin(q1) 0];
P12 = [L2*cos(q2) L2*sin(q2) 0];

% Homogeneous transformation matrices
A01 = [R01(1) R01(4) R01(7) P01(1); R01(2) R01(5) R01(8) P01(2); R01(3) R01(6) R01(9) P01(3); 0 0 0 1]
A12 = [R12(1) R12(4) R12(7) P12(1); R12(2) R12(5) R12(8) P12(2); R12(3) R12(6) R12(9) P12(3); 0 0 0 1]

% Center of mass
MMP01 = simplify(A01*[-L1+Lc1;0;0;1])
MMP02 = simplify(A01*A12*[-L2+Lc2;0;0;1])

% Jacobians and Angular & Translational Velocities for links

syms pL1 pL2 p0 p1 JL1P JL1P1 JL1P2 JL1O JL1O1 JL1O2 JL2P JL2P1 JL2P2 JL2O JL2O1 JL2O2 w1 w2 z dq1 dq2;

z = [0; 0; 1];

JL1P1 = cross(z,[MMP01(1); MMP01(2); MMP01(3)]);
JL1P = [JL1P1(1) 0; JL1P1(2) 0; JL1P1(3) 0]
JL1O1 = z;
JL1O = [z(1) 0; z(2) 0; z(3) 0]

JL2P1 = cross(z,[MMP02(1); MMP02(2); MMP02(3)]);
JL2P = [JL2P1(1) 0; JL2P1(2) 0; JL2P1(3) 0]
JL2O1 = z;
JL2O = [z(1) 0; z(2) 0; z(3) 0]

% Angular velocity

w1 = JL1O * dq1
w2 = JL2O * (dq1 + dq2)

% Translational velocity

dp1 = JL1P * dq1
dp2 = JL2P * dq2

% Inertia tensor
syms IL1 IL2 I0L1 IL1L2 I0L2 R01T R12T R02 R02T m L1 L2
IL1 = [(1/12)*m*L1 0 0; 0 0 0; 0 0 (1/12)*m*L2]
IL2 = [(1/12)*m*L2 0 0; 0 0 0; 0 0 (1/12)*m*L2]

R01T = transpose(R01)
R12T = transpose(R12)
R02 = [cos(q1 + q2) -sin(q1 + q2) 0; sin(q1 + q2) cos(q1 + q2) 0; 0 0 1]
R02T = transpose(R02)

I0L1 = R01T*IL1*R01
I0L1 = simplify(I0L1)
I0L2 = R02T*IL2*R02
I0L2 = simplify(I0L1 + IL1L2)

% Kinetisk energi

syms Ekin(q) ml1 ml2 dp1T dp2T w1T w2T;

dp1T = transpose(dp1)
dp2T = transpose(dp2)
w1T = transpose(w1)
w2T = transpose(w2)

Ekin(q) = 1/2*ml1*dp1T*dp1+1/2*w1T*I0L1*w1+1/2*ml2*dp2T*dp2+1/2*w2T*I0L2*w2
Ekin(q) = simplify(Ekin(q))

pretty(Ekin(q))
latex(Ekin(q))
