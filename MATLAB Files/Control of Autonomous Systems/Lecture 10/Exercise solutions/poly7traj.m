function AC = poly7traj(IC,FC,T)

%POLY7TRAJ Computes coefficients of polynomial for trajectory planning.
%   AC = POLY7TRAJ(IC,FC,T) finds the vector AC of the coefficients of a 
%   polynomial h(t) of degree 7 (h(t)=AC(1)+AC(2)*t+AC(3)*t^2+...+AC(8)*T^7), 
%   with initial conditions being specified in the vector IC=[h(0);h^(1)(0);h^(2)(0);h^(3)(0)], 
%   and final conditions at time t=T specified in the vector 
%   FC=[h(T);h^(1)(T);h^(2)(T);h^(3)(T)].

%   Jerome Jouffroy, January 2023

AC = zeros(8,1);

AC(1) = IC(1); % alpha 0
AC(2) = IC(2); % alpha 1
AC(3) = IC(3)/2; % alpha 2
AC(4) = IC(4)/6; % alpha 3

Tmatrix = [ T^4      ,    T^5  , T^6     ,    T^7    ;
            4*T^4    ,   5*T^5 ,  6*T^6  , 7*T^7     ;
            3*4*T^4  ,4*5*T^5  ,5*6*T^6  ,6*7*T^7    ;
            2*3*4*T^4,3*4*5*T^5,4*5*6*T^6,5*6*7*T^7 ];
        
AC(5:8,1) = Tmatrix\[ FC(1)-(AC(1)+AC(2)*T+AC(3)*T^2+AC(3)*T^3) ;
                            T*(FC(2)-(AC(2)+2*AC(3)*T+3*AC(4)*T^2));
                            T^2*(FC(3)-(2*AC(3)+2*3*AC(4)*T)) ;
                            T^3*(FC(4)-2*3*AC(3))];