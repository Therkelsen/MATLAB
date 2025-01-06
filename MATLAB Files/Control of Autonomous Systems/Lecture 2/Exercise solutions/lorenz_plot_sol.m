
% Plot evolution of state ('trace') for the Lorenz system
% Jerome Jouffroy, September 2024

x_trace = squeeze(out.x); % transform 3D array into 2D array

close all % close all windows

figure, plot3(x_trace(1,:),x_trace(2,:),x_trace(3,:))
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Evolution of state for Lorenz system')
