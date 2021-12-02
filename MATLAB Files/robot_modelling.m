s = settings;
s.matlab.fonts.editor.normal.Size.PersonalValue = 20;
s.matlab.fonts.editor.normal.Name.PersonalValue = 'Calibri';

% Rotation, Translation, Transformation & Center of mass for links
clc
clear
syms MMP01 MMP02 A01 A12 R01 R12 P01 P12 q1 q2 l1 l2 lc1 lc2;

% Rotation matrices
R01 = [cos(q1) -sin(q1) 0; sin(q1) cos(q1) 0; 0 0 1];
R12 = [cos(q2) -sin(q2) 0; sin(q2) cos(q2) 0; 0 0 1];

% Translation vectors
P01 = [l1*cos(q1) l1*sin(q1) 0];
P12 = [l2*cos(q2) l2*sin(q2) 0];

% Homogeneous transformation matrices
A01 = [R01(1) R01(4) R01(7) P01(1); R01(2) R01(5) R01(8) P01(2); R01(3) R01(6) R01(9) P01(3); 0 0 0 1]
A12 = [R12(1) R12(4) R12(7) P12(1); R12(2) R12(5) R12(8) P12(2); R12(3) R12(6) R12(9) P12(3); 0 0 0 1]

% Center of mass
MMP01 = simplify(A01*[-l1+lc1;0;0;1])
MMP02 = simplify(A01*A12*[-l2+lc2;0;0;1])

% Jacobians and Angular & Translational Velocities for links

syms pl1 pl2 p0 p1 Jl1P Jl1P1 Jl1P2 Jl1O Jl1O1 Jl1O2 Jl2P Jl2P1 Jl2P2 Jl2O Jl2O1 Jl2O2 w1 w2 z dq1 dq2;

z = [0; 0; 1];

Jl1P1 = cross(z,[MMP01(1); MMP01(2); MMP01(3)]);
Jl1P = [Jl1P1(1) 0; Jl1P1(2) 0; Jl1P1(3) 0]
Jl1O1 = z;
Jl1O = [z(1) 0; z(2) 0; z(3) 0]

Jl2P1 = cross(z,[MMP02(1); MMP02(2); MMP02(3)]);
Jl2P = [Jl2P1(1) 0; Jl2P1(2) 0; Jl2P1(3) 0]
Jl2O1 = z;
Jl2O = [z(1) 0; z(2) 0; z(3) 0]

% Angular velocity
w1 = Jl1O * dq1
w2 = Jl2O * (dq1 + dq2)

% Translational velocity
dp1 = Jl1P * dq1
dp2 = Jl2P * dq2

% Inertia tensor
syms I Il1 Il2 I0l1 Il1l2 I0l2 R01T R12T R02 R02T m1 m2 l1 l2
%Il1 = [(1/12)*m1*l1 0 0; 0 0 0; 0 0 (1/12)*m1*l2]
%Il2 = [(1/12)*m2*l2 0 0; 0 0 0; 0 0 (1/12)*m2*l2]

R01T = transpose(R01)
R12T = transpose(R12)
R02 = [cos(q1 + q2) -sin(q1 + q2) 0; sin(q1 + q2) cos(q1 + q2) 0; 0 0 1]
R02T = transpose(R02)

I0l1 = R01T*I*R01
I0l1 = simplify(I0l1)
I0l2 = R02T*I*R02
I0l2 = simplify(I0l2)

% Kinetisk energi

syms Ekin(q) ml1 ml2 dp1T dp2T w1T w2T;

dp1T = transpose(dp1)
dp2T = transpose(dp2)
w1T = transpose(w1)
w2T = transpose(w2)

Ekin(q) = (1/2*ml1*dp1T*dp1)+(1/2*w1T*I0l1*w1)+(1/2*ml2*dp2T*dp2)+(1/2*w2T*I0l2*w2)

pretty(Ekin(q))
