function AC = poly5traj(IC,FC,T)

%POLY5TRAJ Computes coefficients of polynomial for trajectory planning.
%   AC = POLY5TRAJ(IC,FC,T) finds the vector AC of the coefficients of a 
%   polynomial h(t) of degree 5 (h(t)=AC(1)+AC(2)*t+AC(3)*t^2+...+AC(6)*T^5), 
%   with initial conditions being specified in the vector IC=[h(0);h^(1)(0);h^(2)(0)], 
%   and final conditions at time t=T specified in the vector 
%   FC=[h(T);h^(1)(T);h^(2)(T)].

%   Jerome Jouffroy, January 2023

AC = zeros(6,1);

AC(1) = IC(1); % alpha 0
AC(2) = IC(2); % alpha 1
AC(3) = IC(3)/2; % alpha 2

Tmatrix = [ T^3      ,    T^4  , T^5        ;
            3*T^3    ,   4*T^4 ,  5*T^5     ;
            2*3*T^3  ,3*4*T^4  ,4*5*T^5   ];
        
AC(4:6,1) = inv(Tmatrix)*[ FC(1)-(AC(1)+AC(2)*T+AC(3)*T^2) ;
                            T*(FC(2)-(AC(2)+2*AC(3)*T));
                            T^2*(FC(3)-2*AC(3)) ];