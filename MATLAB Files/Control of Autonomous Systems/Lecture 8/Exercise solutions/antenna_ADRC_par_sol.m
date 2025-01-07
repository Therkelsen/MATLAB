% parameter file for the ADRC of the satellite antenna
% Jerome Jouffroy, October 2023

% system parameters
a_m = 1.71;
Kl = 100;
Kg = 0.1;

% state-space representation
A = [ 0 1 ; 0 -a_m ];
B = [ 0 ; Kl*Kg ];
C = [ 1 , 0 ];
D = 0;

% ADRC parameters
b = Kl*Kg;
A_cont = [ 0 , 1 ; 0 , 0 ];
B_cont = [ 0 ; 1];
pole_cont = -1;
p_cont = pole_cont*ones(2,1);
K = acker(A_cont,B_cont,p_cont);
k1 = K(1);
k2 = K(2);

% ESO tuning
p_ESO = 10*pole_cont*ones(3,1);
A_ESO = [ 0 , 1 , 0 ; 0 , 0 , 1 ; 0 , 0 , 0 ];
B_ESO = [ 0; b ; 0 ];
C_ESO = [ 1 , 0 , 0 ];
L_ESO = acker(A_ESO',C_ESO',p_ESO)';