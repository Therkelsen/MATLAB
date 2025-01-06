clc; clear; close all;
format compact;

% Hot air balloon parameters
% Jerome Jouffroy, November 2024

% state-space representation
A = [ -1 0 0 ; 0 0 1 ; 1 0 -1/2 ];
B = [ 1 ; 0 ; 0 ];
C = [ 0 , 1 , 0 ];

% computations towards obtaining a CCF after change of coordinates.
Wc = [ B , A*B , A^2*B ];
det_Wc = det(Wc)

Bt = [ 0 ; 0 ; 1 ];

T1 = Bt'*inv(Wc)
T2 = T1*A;
T3 = T2*A;

Tcoord = [ T1 ; T2 ; T3 ] % matrix T for change of coordinates

At = Tcoord*A*inv(Tcoord);

% selection of the a_i's parameters from matrix Atilde
a_vect = -At(3,:)';

%------------------
% motion planning.
%------------------
T = 600; % horizon time
% intial conditions
h0 = 0;
v0 = 0;
a0 = 0;

% final conditions
hT = 100;
vT = 0;
aT = 0;

% calculations polynomial coefficients
alpha_vect = poly5traj([h0;v0;a0],[hT;vT;aT],T);