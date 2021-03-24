syms A V L a theta t1 t2 t3 t4 L1 L2 L3 result;
theta =             [t1 t2 t3 t4];
L =                 [L1 L2 L3];
T01 =                 [cos(theta(1)) -sin(theta(1)) 0 0; sin(theta(1)) cos(theta(1)) 0 0; 0 0 1 0; 0 0 0 1];
T12 =                 [cos(theta(2)) -sin(theta(2)) 0 L(1); sin(theta(2)) cos(theta(2)) 0 0; 0 0 1 0; 0 0 0 1];
T23 =                 [cos(theta(3)) -sin(theta(3)) 0 L(2); sin(theta(3)) cos(theta(3)) 0 0; 0 0 1 0; 0 0 0 1];

T01*T12*T23
