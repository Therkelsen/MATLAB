clc; clear;close all

%% Dynamic model of hexacopter

% PID GAINS
% Altitude: 1,0,8
% Roll : 1,0,5
% Pitch : 1,0,5
% Yaw : 1,0,3

% Simulation parameters:
    g = 9.82; % Acceleration of gravity [m/s^2]
    % Aero parameters
    K_phi = 0.005; % Roll damping [N*m*s/rad]
    K_theta = 0.005; % pitch damping [N*m*s/rad]
    K_psi = 0.002; % Yaw damping [N*m*s/rad]
    K_x = 0.3; % Forwards damping [N*s/m]
    K_y = 0.3; % Sideways damping [N*s/m]
    K_z = 0.2; % Vertical damping [N*s/m]
    
    Simulation_params = [g,K_phi,K_theta,K_psi,K_x,K_y,K_z];

% Body parameters:
    L = 0.275; % Arm length [m]
    b = 1.177e-5; % Thrust coefficient [N/m^2]
    d = 1.855e-7; % Drag coefficient
    m = 2.5; % Mass [kg]
    Ix = 100*0.00915; % moment of inertia x [kg*m^2]
    Iy = 100*0.00915; % moment of inertia y [kg*m^2]
    Iz = 0.01187; % moment of inertia z [kg*m^2]

    % Rotor parameters:
    m_rotor = 0.1; % Mass of propeller + motor [kg]
    r_rotor = 0.1; % radius of propeller [m]
    J_T = 1/2*m_rotor*r_rotor; % Inertia of each rotor [kg*m^2]

    body_params = [L,b,d,m,Ix,Iy,Iz,J_T];

% Thrust matrix:
T = [1, 1, 1, 1, 1, 1;
     -L, L, L/2, -L/2, -L/2, L/2;
     0, 0, L*sqrt(3)/2, -L*sqrt(3)/2, L*sqrt(3)/2, -L*sqrt(3)/2;
     -d/b, d/b, -d/b, -d/b, d/b, d/b];

T_inv = pinv(T);

%% Attitude and altitude controller

% Sliding mode control parameters:
    k = 10; 
    delta = 3;

    control_params = [k,delta];

% Fault detection parameters:
    A = [zeros(4),eye(4);
         zeros(4),zeros(4)];
  K = [eye(4)*1.1;eye(4)*0.3025]*10;
   % C = [eye(4),zeros(4)];
   C = eye(8);
    B = [zeros(4);
         0,-L/Ix,0,L/Ix;  %0,-L/Iy,0,L/Ix?
         -L/Iy,0,L/Iy,0;
         d/(b*Iz),-d/(b*Iz),d/(b*Iz),-d/(b*Iz);
         1/m,1/m,1/m,1/m;];
    desired_poles = [-10+2i,-10-2i,-40+2i,-40-2i,-30+2i,-30-2i,-20+2i,-20-2i];
    K = place(A', C', desired_poles)'
    