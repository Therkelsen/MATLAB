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