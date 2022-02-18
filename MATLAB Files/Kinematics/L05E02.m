syms K
syms v0
syms theta
theta = degtorad(30)
K = [(1/sqrt(3)) -(1/sqrt(3)) (1/sqrt(3))]
v0 = 1-cos(theta)
syms RK
RK = [(K(1)*K(1)*v0+cos(theta)) (K(1)*K(2)*v0-K(3)*sin(theta)) (K(1)*K(3)*v0+K(2)*sin(theta)); (K(1)*K(2)*v0+K(3)*sin(theta)) (K(2)*K(2)*v0+cos(theta)) (K(2)*K(3)*v0-K(1)*sin(theta)); (K(1)*K(3)*v0-K(2)*sin(theta)) (K(2)*K(3)*v0+K(1)*sin(theta)) (K(3)*K(3)*v0+cos(theta))]
