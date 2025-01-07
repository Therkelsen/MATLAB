% plot a simple 2D parabola and finds its minimum using quadprog
% Jerome Jouffroy, October 2024

close all

[X1,X2] = meshgrid(-10:1:10,-10:1:10);
J = (X1-1).^2 + (X2-2).^2;
figure, surf(X1,X2,J)
title('2D parabola')
xlabel('x_1')
ylabel('x_2')
hold on

% find an optimum without constraints on x1 and x2
H = 2*eye(2);
F = [ -2 ; -4 ]';

options = optimoptions('quadprog','Display','off');
xopt = quadprog(H,F)
x1opt = xopt(1,1);
x2opt = xopt(2,1);
Jopt = (x1opt-1)^2 + (x2opt-2)^2;
plot3(x1opt,x2opt,Jopt,'r.','MarkerSize',20)

% find an optimum under constraints outside of global minimum
lb = [ -5 ; -5 ];
ub = [ 5 ; -3 ];
xsubopt = quadprog(H,F,[],[],[],[],lb,ub,[])
x1subopt = xsubopt(1,1);
x2subopt = xsubopt(2,1);
Jsubopt = (x1subopt-1)^2 + (x2subopt-2)^2;
plot3(x1subopt,x2subopt,Jsubopt,'g.','MarkerSize',20)
view(110,50)
hold off

