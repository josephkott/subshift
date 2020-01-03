%% Test MEX scanner
clc; clear

omega = -1;
L1 = pi / 2; L2 = pi / 2;
params = [omega L1 L2];

intervals = 512;

% Scan plane grid
M = 512; N = 512;
u_span = [-4 +4];
du_span = [-6 +6];

% \mathcal{U}_{\pi}^+ set
x_span_plus = [0 +pi];

% \mathcal{U}_{\pi}^- set
x_span_minus = [0 -pi];

% With time measurements
tic

U_plus_set = f_scan(params, x_span_plus, u_span, du_span, M, N, intervals);

t = toc;
figure
subplot(2, 2, 1)
fprintf('Elapsed f_scan: %f\n', t)
plot_scan(u_span, du_span, U_plus_set, false) % not create new figure
title(sprintf('Elapsed f_scan: %f', t), 'Interpreter', 'none')

tic

U_minus_set = f_scan(params, x_span_minus, u_span, du_span, M, N, intervals); 

t = toc;
fprintf('Elapsed f_scan: %f\n', t)
subplot(2, 2, 2)
plot_scan(u_span, du_span, U_minus_set, false) % not create new figure
title(sprintf('Elapsed f_scan: %f', t), 'Interpreter', 'none')

%% MATLAB scanner implementation

u_step = (u_span(2) - u_span(1)) / (M - 1);
du_step = (du_span(2) - du_span(1)) / (N - 1);

u_grid = u_span(1) : u_step : u_span(2);
du_grid = du_span(1) : du_step : du_span(2);

U_plus_set_ml = zeros(N, M);

tic

for i = 1 : length(du_grid)
	for j = 1 : length(u_grid)
		[X, U, is_infinite] = f_solve(params, x_span_plus, [u_grid(j) du_grid(i)], intervals);
		U_plus_set_ml(i, j) = is_infinite;
	end
end

t = toc;
fprintf('Elapsed f_solve + loop: %f\n', t)
subplot(2, 2, 3)
plot_scan(u_span, du_span, U_plus_set_ml, false) % not create new figure
title(sprintf('Elapsed f_scan: %f', t), 'Interpreter', 'none')

U_minus_set_ml = zeros(N, M);

tic

for i = 1 : length(du_grid)
	for j = 1 : length(u_grid)
		[X, U, is_infinite] = f_solve(params, x_span_minus, [u_grid(j) du_grid(i)], intervals);
		U_minus_set_ml(i, j) = is_infinite;
	end
end

t = toc;
fprintf('Elapsed f_solve + loop: %f\n', t) % not create new figure
subplot(2, 2, 4)
plot_scan(u_span, du_span, U_minus_set_ml)
title(sprintf('Elapsed f_scan: %f', t), 'Interpreter', 'none')