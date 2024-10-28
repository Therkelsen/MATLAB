% Define the two functions
f1 = @(N) (log2(N)).^2 .* N.^(3/2);  % O_1: (log2(N))^2 * N^(3/2)
f2 = @(N) N.^2;                      % O_2: N^2

% Define a function for the difference (f1 - f2) to find where it equals 0
diff_fun = @(N) f1(N) - f2(N);

% Set options for fsolve (optional)
options = optimset('Display','off', 'TolX', 1e-8);  % Increase accuracy

% Set a larger initial guess near the expected intersection
N_guess_large = 65000;  % Initial guess near the expected intersection

% Use fsolve to find the intersection
N_intersect_large = fsolve(diff_fun, N_guess_large, options);

% Display the result
disp(['Intersection at N = ', num2str(N_intersect_large)]);

% Plot the functions to visualize
N_values = linspace(1, 100000, 1000);  % Adjust range as needed
figure;
plot(N_values, f1(N_values), 'r-', 'LineWidth', 2); hold on;
plot(N_values, f2(N_values), 'b-', 'LineWidth', 2);
plot(N_intersect_large, f1(N_intersect_large), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g');  % Mark the intersection
xlabel('N');
ylabel('Function Value');
legend('O_1(N)', 'O_2(N)', 'Intersection');
title('Intersection of O_1 and O_2');
grid on;
