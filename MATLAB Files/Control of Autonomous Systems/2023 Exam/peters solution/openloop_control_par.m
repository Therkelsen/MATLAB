clc, clear

% Constants a = [q1, q2, q3, J]

q1 = 1;
q2 = 10;
q3 = 0.1;
J = 40;

a = [q1, q2, q3, J];

x0 = [2; 200];

% x2_star = x0(2);
% x1_star = q3/q2 * x2_star;
% u_star = q1 * x1_star * x2_star;

% state-space representation
% We rewrite the cancelling controller
% Hence the v replaces the u.
A = [0, 0; q2/J, -q3/J];
B = [ 1 ; 0];



% computations towards obtaining a CCF after change of coordinates.
Wc = [ B , A*B];


Bt = [ 0 ; 1 ];

%Since number og states is 2
T1 = Bt'*inv(Wc)
T2 = T1*A;

Tcoord = [ T1 ; T2] % matrix T for change of coordinates

At = Tcoord*A*inv(Tcoord);

a_vect = -At(2,:)';

p0 = 200;
v0 = 0;

pT = 300;
vT = 0;

% motion planning
H = 10;

% alpha_vect = poly3traj([p0;v0],[pT;vT],H);

x0 = [ 2 ; 200 ];
xT = [ 3 ; 300 ];
H = 10;

z0 = Tcoord*x0;
zT = Tcoord*xT;

alpha_vect = poly3traj(z0,zT,H);


function AC = poly3traj(IC,FC,T)

%POLY3TRAJ Computes coefficients of polynomial for trajectory planning.
%   AC = POLY3TRAJ(IC,FC,T) finds the vector AC of the coefficients of a 
%   polynomial h(t) of degree 3 (h(t)=AC(1)+AC(2)*t+AC(3)*t^2+AC(4)*t^3), 
%   with initial conditions being specified in the vector IC=[h(0);h_dot(0)], 
%   and final conditions at time t=T specified in the vector 
%   FC=[h(T);h_dot(T)].

%   Jerome Jouffroy, October 2009

AC = zeros(4,1);

AC(1) = IC(1);
AC(2) = IC(2);

Tmatrix = [ T^2    ,   T^3  ;
            2*T^2  , 3*T^3 ]/T^2;

FCandAC = [FC(1)-(AC(1)+AC(2)*T);T*(FC(2)-AC(2))];

AC(3:4,1) = inv(Tmatrix)/T^2*FCandAC;
end
