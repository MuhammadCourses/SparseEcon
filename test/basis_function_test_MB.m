addpath(genpath('/Users/muhammadbashir/GitHub/MuhammadCourses/SparseEcon'))

% Define parameters
n = 10; % Number of basis functions
nodes = linspace(0, 1, 100); % Points at which to evaluate the basis functions
region_min = 0; % Minimum of the region
region_max = 1; % Maximum of the region
type = 'cheb'; % Type of basis functions
X = linspace(0, 1, 2*n); % Vector of nodal basis points

% Call basis_function to get basis functions T
[T, T_prime] = basis_function(n, nodes, region_min, region_max, type, X);

% Define the function f(x)
f = @(x) (4/5) * (sin(pi * x) + (1/2) * sin(2 * pi * x));

% Evaluate f(x) at the nodes
f_values = f(nodes);

% Approximate f(x) using the basis functions
% Solve for coefficients a such that T' * a = f_values
a = T' \ f_values';

% Reconstruct the approximation of f(x) using the basis functions
f_approx = T' * a;

% Plot the original function and the approximation
figure;
hold on;
plot(nodes, f_values, 'r', 'DisplayName', 'Original Function f(x)');
plot(nodes, f_approx, 'b--', 'DisplayName', 'Approximation using Basis Functions');
hold off;
xlabel('Nodes');
ylabel('Function Value');
title('Function Approximation using Nodal Basis Functions');
legend show;

% Display the basis functions
disp('Basis Functions (T):');
disp(T);


% here let us practice use of basis_fun_irf to get the coefficients and projections of the function
% Define parameters
H = 5; % Number of basis functions
N = 2; % Number of functions
t = linspace(0, 1, 100); % Time points or nodes

% Define the functions Y
Y = zeros(length(t), N);
Y(:, 1) = sin(pi * t); % First function
Y(:, 2) = cos(pi * t); % Second function

% Get coefficients
coeffs = basis_fun_irf(Y, [], H, N, 'cheb', t, 'get_coefficient');
% Reconstruct functions
reconstructed_funcs = basis_fun_irf([], coeffs, H, N, 'cheb', t, 'get_function');


% Plot the original and reconstructed functions
figure;
subplot(2, 1, 1);
plot(t, Y(:, 1), 'r', 'DisplayName', 'Original sin(pi * t)');
hold on;
plot(t, reconstructed_funcs(:, 1), 'b--', 'DisplayName', 'Reconstructed sin(pi * t)');
hold off;
xlabel('t');
ylabel('Function Value');
title('Original and Reconstructed Functions');
legend show;

subplot(2, 1, 2);
plot(t, Y(:, 2), 'r', 'DisplayName', 'Original cos(pi * t)');
hold on;
plot(t, reconstructed_funcs(:, 2), 'b--', 'DisplayName', 'Reconstructed cos(pi * t)');
hold off;
xlabel('t');
ylabel('Function Value');
legend show;


