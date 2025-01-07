% DT double-integrator MPC example
% Jerome Jouffroy, October 2023

A = [ 1 1 ; 0 1 ];
B = [ 0 ; 1 ];
C = [ 1 0 ];
x0 = [ 10 ; 0 ];

Q = C'*C;
R = 1/10;

Qb = blkdiag(Q,Q,Q);
Rb = blkdiag(R,R,R);

Umin = -1*ones(3,1);
Umax = 1*ones(3,1);

Ab = [ eye(2) ; A ; A^2 ];
Bb = [ zeros(2,3) ; B zeros(2,2) ; A*B B zeros(2,1) ];

H = 2*(Bb'*Qb*Bb + Rb);
preF = Ab'*Qb*Bb;

N = 41; % number of iterations

xk = x0;

x_trace = zeros(2,N);
u_trace = zeros(1,N);

options = optimoptions('quadprog','Display','off');

for k=0:N-1
    % compute F matrix
    %F = 2*(xk'*Ab'*Q4*Bb)';
    F = 2*(xk'*preF)'; % more efficient calculation
    % compute optimal vector Ustar
    Ustar = quadprog(H,F,[],[],[],[],Umin,Umax,[],options);
    % select current input
    uk = Ustar(1);
    % apply control input
    xkp1 = A*xk + B*uk;
    % store value of state and control input in trace
    x_trace(:,k+1) = xk;
    u_trace(1,k+1) = uk;
    % prepare next iteration
    xk = xkp1;
end

close all

figure, plot(0:N-1,x_trace(1,:),0:N-1,x_trace(2,:),'LineWidth',2)
grid on
title('Evolution of state')
xlabel('iterations k')

figure, plot(0:N-1,u_trace,'LineWidth',2)
grid on
title('Evolution of control input')
xlabel('iterations k')
    