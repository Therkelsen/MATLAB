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

Jl1P1 = cross(z,[MMP01(1)-0; MMP01(2)-0; MMP01(3)-0]);
Jl1P = [Jl1P1(1) 0; Jl1P1(2) 0; Jl1P1(3) 0]
Jl1O1 = z;
Jl1O = [z(1) 0; z(2) 0; z(3) 0]

Jl2P1 = cross(z,[MMP02(1)-0; MMP02(2)-0; MMP02(3)-0]);
Jl2P2 = cross(z,[MMP02(1)-P01(1); MMP02(2)-P01(2); MMP02(3)-P01(3)]);
Jl2P = [Jl2P1(1) Jl2P2(1); Jl2P1(2) Jl2P2(2); Jl2P1(3) Jl2P2(3)]
Jl2O1 = z;
Jl2O = [Jl2O1(1) Jl1O1(1); Jl2O1(2) Jl1O1(2); Jl2O1(3) Jl1O1(3)]

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


% Stedvektor til massemidtpunkt
p0c1 = [lc1*cos(q1); lc1*sin(q1); 0]
p0c2 = [lc2*cos(q1 + q2) + l1*cos(q1); lc2*sin(q1 + q2) + l1*sin(q1); 0]


% Potentiel energi
syms g ml1 ml2 
g0 = [0; -g; 0]
g0T = transpose(g0)

Epot = -(ml1 * g0T *p0c1 + ml2 * g0T * p0c2)


% Kinetisk energi
syms Ekin(q) dp1T dp2T w1T w2T dq I1l1 I2l2;
dq = [dq1; dq2]
dqT = transpose(dq)
dp1T = transpose(dp1)
dp2T = transpose(dp2)
w1T = transpose(w1)
w2T = transpose(w2)
Jl1PT = transpose(Jl1P)
Jl1OT = transpose(Jl1O)
Jl2PT = transpose(Jl2P)
Jl2OT = transpose(Jl2O)

%Ekin(q) = (1/2*ml1*dp1T*dp1)+(1/2*w1T*I0l1*w1)+(1/2*ml2*dp2T*dp2)+(1/2*w2T*I0l2*w2)

Ekinl1(q) = 1/2 * ml1 * (dqT * Jl1PT) * (Jl1P * dq) + 1/2 * (dqT * Jl1OT) * R01 * I1l1 * R01T * (Jl1O * dq);
Ekinl2(q) = 1/2 * ml2 * (dqT * Jl2PT) * (Jl2P * dq) + 1/2 * (dqT * Jl2OT) * R02 * I2l2 * R02T * (Jl2O * dq);
Ekin(q) = Ekinl1(q) + Ekinl2(q)
% Følgende ligninger viser ligningerne ovenfor er rigtige. Tjek dem mod lektion 11 slide 33.
%Ekinl1Test(q) = ml1 * (Jl1PT) * (Jl1P) +(Jl1OT) * R01 * I1l1 * R01T * (Jl1O)
%Ekinl2Test(q) = ml2 * (Jl2PT) * (Jl2P) + (Jl2OT) * R02 * I2l2 * R02T * (Jl2O)

pretty(Ekin(q))


%% 
% Trans-hastigheder
syms l1 s1 c1 c2 lc1 lc2 c12 s12 voc2 qp qp1 qp2 nul A B C res voc1;

qp = [qp1; qp2]
nul = [0; 0; 1]
A = [l1*c1 + lc2*c12; l1*s1 + lc2*s12; 0]
B = [l1*c1; l1*s1; 0]
voc2 = [cross(nul,A) cross(nul,(A-B))]
res = voc2*qp

C = [-l1*s1 0; lc1*c1 0; 0 0]
voc1 = C*qp


%%
% Bevægelsesligninger
% Ved brug af Euler-Lagrance modellering
syms ddq1 ddq2 tau1 tau2

L = Ekin - Epot
% part1 skal manuelt differentieres i forhold til tid efter partial differention
% part1 = diff(L,dq)
part11 = diff(L,dq1)
% simplify(part11) = I1l1*dq1 + I2l2*dq1 + I2l2*dq2 + dq1*l1^2*ml2 + dq1*lc1^2*ml1 + dq1*lc2^2*ml2 + dq2*lc2^2*ml2 + 2*dq1*l1*lc2*ml2*cos(q2) + dq2*l1*lc2*ml2*cos(q2)
% manuel simplify = dq1*(I1l1 + I2l2 + lc1^2*ml1 + ml2*(l1^2 + lc2^2 + 2*l1*lc2*cos(q2))) + dq2*(I2l2 + ml2*(lc2^2 + l1*lc2*cos(q2)))
part12 = diff(L,dq2)
% simplify(part12) = I2l2*dq1 + I2l2*dq2 + dq1*lc2^2*ml2 + dq2*lc2^2*ml2 + dq1*l1*lc2*ml2*cos(q2)
% manual simplify = dq1*(I2l2 + lc2^2*ml2 + l1*lc2*ml2*cos(q2)) + dq2*(I2l2 + lc2^2*ml2)
% part2 = diff(L,q)
part21 = diff(L,q1)
part22 = diff(L,q2)

% Lav manuel diff
% simplify(part11) = I1l1*dq1 + I2l2*dq1 + I2l2*dq2 + dq1*l1^2*ml2 + dq1*lc1^2*ml1 + dq1*lc2^2*ml2 + dq2*lc2^2*ml2 + 2*dq1*l1*lc2*ml2*cos(q2) + dq2*l1*lc2*ml2*cos(q2)
part11d = I1l1*ddq1 + I2l2*ddq1 + I2l2*ddq2 + ddq1*l1^2*ml2 + ddq1*lc1^2*ml1 + ddq1*lc2^2*ml2 + ddq2*lc2^2*ml2 + 2*ddq1*l1*lc2*ml2*(-sin(q2)) + ddq2*l1*lc2*ml2*(-sin(q2))
%simplify(part12) = I2l2*dq1 + I2l2*dq2 + dq1*lc2^2*ml2 + dq2*lc2^2*ml2 + dq1*l1*lc2*ml2*cos(q2)
part12d = I2l2*ddq1 + I2l2*ddq2 + ddq1*lc2^2*ml2 + ddq2*lc2^2*ml2 + ddq1*l1*lc2*ml2*(-sin(q2))

full1 = part11d - part21 == tau1
full2 = part12d - part22 == tau2

Complete = [part11d; part12d] - [part21; part22] == [tau1;tau2]

ddq1 = isolate(full1, ddq1)
ddq2 = isolate(full2, ddq2)
% Disse burde virke efter at part11d og part12d er lavet ordentligt (læs
% relevant OBS ved disse)
%dq1 = isolate(full1, dq1) 
%dq2 = isolate(full2, dq2)


% % Eksempel 1
% syms L m k x dx t
% L = 1/2*m*dx^2 - 1/2*k*x^2
% part1 = diff(L,dx)
% part2 = diff(L,x)
% full = part1 - part2
% 
% 
% % Eksempel 2
% syms K1 K2 o1 o2 I1 I2 do1 do2
% kin = 1/2 * I1 * do1^2 + 1/2 * I2 * do2^2
% pot = 1/2 * K1 * o1^2 + 1/2 * K2* (o1 - o2)^2
% L = kin-pot
% 
% part11 = diff(L,do1)
% part21 = diff(L,o1)
% full1 = part11 - part21
% 
% part12 = diff(L,do2)
% part22 = diff(L,o2)
% full2 = part12 - part22

% clear
% % Ekspempel 3
% syms M m dxc ddxc xc m l g do ddo o Jp d dt
% L = 1/2*(M+m)*dxc^2 + m*l*dxc*do*cos(o) + 1/2*(m*l^2+Jp)*do^2 + m*g*l*cos(o)
% part11 = diff(L,dxc)
% part12 = diff(L,do)
% part21 = diff(L,xc)
% part22 = diff(L,o)
% % Først skal part11 og part12 differentieres i forhold til tiden (måske
% % bare sæt d/dt foran symbols. Det har Christoffer gjort i slidene
% Complete = [part11;part12] - [part21; part22]
