% Definition of parameters for the identification of the Piper model
% Jerome Jouffroy, September 2024

% definition of range for different omegas
omega = 1:9;
omega_low = 0.1*omega;
omega_mid = omega;
omega_high = 10*omega;

omegas = [ omega_low , omega_mid , omega_high ];

% definition of randomized phases
phis = 2*pi*rand(1,length(omegas)); 