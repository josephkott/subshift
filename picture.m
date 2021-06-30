%% tmp
clc
clear

L1 = 2;
L2 = 1;
L = L1 + L2;

params = [-1 L1 L2];
u_span = [-2 +2]; du_span = [-4 +4];
M = 256; N = 256;

U_plus_set = f_scan(params, [0, L1], u_span, du_span, M, N, 1024);
U_minus_set = f_scan(params, [0 -L], u_span, du_span, M, N, 1024);
plot_scan(u_span, du_span, U_plus_set  + U_minus_set, true, 'gray')

%% Parameters
clc
clear

L1 = 2;
L2 = 1;
L = L1 + L2;

params = [-1 L1 L2];
u_span = [-2 +2]; du_span = [-4 +4];
M = 1024; N = 1024;

%% Picture blocks

fprintf('U_plus_set, P1U_plus_set\n')
U_plus_set = f_scan(params, [0, L1], u_span, du_span, M, N, 1024);
P1U_plus_set = U_plus_set(end:-1:1, :);
plot_scan(u_span, du_span, U_plus_set, true, 'gray')
plot_scan(u_span, du_span, P1U_plus_set, true, 'gray')

fprintf('U_minus_set\n')
U_minus_set = f_scan(params, [0 -L], u_span, du_span, M, N, 1024);
plot_scan(u_span, du_span, U_minus_set, true, 'gray')

fprintf('P1Islands_set\n')
U_half_minus_minus_set = f_scan(params, [L1 -L], u_span, du_span, M, N, 1024);
U_half_minus_set_with_bounds = f_scan_with_bounds(params, [L1 0], u_span, du_span, ...
	                                              M, N, max(u_span), max(du_span), 1024);
P1Islands_set = U_half_minus_minus_set + U_half_minus_set_with_bounds ~= 0;
plot_scan(u_span, du_span, P1Islands_set, true, 'gray')

fprintf('PIslands_set\n')
U_minus_minus_set = f_scan(params, [L -L], u_span, du_span, M, N, 1024);
U_minus_set_with_bounds = f_scan_with_bounds(params, [L 0], u_span, du_span, ...
											 M, N, max(u_span), max(du_span), 1024);
PIslands_set = U_minus_minus_set + U_minus_set_with_bounds ~= 0;
plot_scan(u_span, du_span, PIslands_set, true, 'gray')

%% Separatrix
u = u_span(1):0.1:u_span(2);
du = sqrt(0.5 * (u .^ 4) - params(1) * (u .^ 2));
figure('Position', [100, 100, 350, 300])
plot(u, du, u, -du, 'LineWidth', 2, 'Color', 'red')
axis([u_span(1) u_span(2) du_span(1) du_span(2)])

%%
M = 1024; N = 1024;
U_set = f_scan_with_energy_bound(params, [0 (2 * L)], u_span, du_span, M, N, 50, 1024);
% U_set = f_scan(params, [0 (2 * L)], u_span, du_span, M, N, 1024);
plot_scan(u_span, du_span, U_set, true)

%% Islands map
P1_islands = map_islands(params, u_span, du_span, M, N, U_plus_set, U_minus_set, L1);
plot_scan(u_span, du_span, P1U_plus_set + P1_islands, true)

%%
U_minus_minus_set = f_scan(params, [L -L], u_span, du_span, M, N, 1024);
plot_scan(u_span, du_span, U_minus_minus_set, true)

%%
U_minus_set_with_bounds = f_scan_with_bounds(params, [L 0], u_span, du_span, ...
											 M, N, max(u_span), max(du_span), 1024);
PIslands_set = U_minus_minus_set + U_minus_set_with_bounds ~= 0;
plot_scan(u_span, du_span, PIslands_set, true)

%%
clc
clear

L1 = 2;
L2 = 0.1;
L = L1 + L2;

params = [-1 L1 L2];
u_span = [9.97 +9.99]; du_span = [-71.35 -71];
M = 128; N = 128;

U_plus_set = f_scan(params, [0, L1], u_span, du_span, M, N, 1024);
U_minus_set = f_scan(params, [0 -L], u_span, du_span, M, N, 1024);
plot_scan(u_span, du_span, U_plus_set + U_minus_set, true)