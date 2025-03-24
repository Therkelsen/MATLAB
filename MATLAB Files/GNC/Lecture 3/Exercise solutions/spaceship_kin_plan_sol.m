% motion plannning for theta angle from 0 to 360 degrees in 10 seconds
% Jerome Jouffroy, January 2024


%------------------
% motion planning.
%------------------
T = 10; % horizon time
% intial conditions
theta0 = 0;
theta_dot0 = 0;
theta_ddot0 = 0;

% final conditions
thetaT = 3*360;
theta_dotT = 0;
theta_ddotT = 0;

% calculations polynomial coefficients
alpha_vect = poly5traj([theta0;theta_dot0;theta_ddot0],[thetaT;theta_dotT;theta_ddotT],T);

