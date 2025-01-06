% parameters and motion planning for the MSD system
% Jerome Jouffroy, November 2024

% parameters
m = 10;
k = 5;
d = 0.1;

a_vect = [ k ; d ; m ];

% setting up initial and final conditions
p0 = 0;
v0 = 0;
a0 = 0;
j0 = 0; % jerk: derivative of acceleration

pT = 30;
vT = 0;
aT = 0;
jT = 0;

% motion planning
H = 10;

% motion planning with different orders: one can see different levels of
% smoothness of control input ud(t) according which order is used
% order 3
alpha_vect = poly3traj([p0;v0],[pT;vT],H);
% order 5
%alpha_vect = poly5traj([p0;v0;a0],[pT;vT;aT],H);
% order 7
%alpha_vect = poly7traj([p0;v0;a0;j0],[pT;vT;aT;jT],H);

ustar = k*pT;