% 0.1
clc; clear; close all;
syms u

y = 3*sin(u)

disp('y for different input values u.')
for c = 1:10
    result = double(subs(y, u, c));
    % fprintf('y(%d) = 3 * sin(%d) -> %.4f\n', c, c, result);  % Print formatted output
end


% 0.2

disp('Values of y for an interval of u ∈ [0, 2π].')
for c = 1:0.1:2*pi
    result = double(subs(y, u, c));
    % fprintf('y(%.2f) = 3 * sin(%.2f) -> %.4f\n', c, c, result);  % Print formatted output
end

% Plot y(u) for u from 0 to 2*pi
u_plot = 0:0.01:2*pi;  % Finer resolution for the plot
y_plot = double(subs(y, u, u_plot));

% figure;
% plot(u_plot, y_plot, 'LineWidth', 2);
% xlabel('u');
% ylabel('y(u) = 3 * sin(u)');
% title('Plot of y(u) = 3 * sin(u)');
% grid on;

% 0.3

fprintf('rotated vector w')
v = transpose([1, 2])

syms theta
R = [cos(theta), -sin(theta);
     sin(theta), cos(theta)]

R = subs(R, theta, 2*pi/3)

w = R*v

% Exercise 1

m = 1;
k = 2;
b = 2;

% Exercise 2

