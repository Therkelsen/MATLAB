% estimation of the transder function of the Piper example
% Jerome Jouffroy, September 2024

% definition of the number of poles and zeros
np = 4;
nz = 2;
% np = 3;
% nz = 1; % what if we make a mistake in the number of poles and zeros?

% group input and output logs into 'data'
Ts = 0.01;
data = iddata(out.y,out.u,Ts);

% estimation of the transfer function P_hat (ie the estimate of the actual
% P)
P_hat = tfest(data,np,nz)

figure, bode(P,'r',P_hat,'b--')
title('comparison between P and its estimate')
grid on