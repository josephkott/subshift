%% Test \mathcal{U}_{L_1}^+ borders
clc; clear

intervals = 512;
M = 512; N = 512;
u_span = [-5 +5];
du_span = [-8 +8];

%% Test #1
omega = -1;
L1 = pi / 2; L2 = pi / 2;
params = [omega L1 L2];

x_span_plus = [0, L1 + L2];
U_plus_set = f_scan(params, x_span_plus, u_span, du_span, M, N, intervals);

figure('Position', [100, 100, 350, 300])
axis([u_span du_span])
hold on

plot_scan(u_span, du_span, U_plus_set, false)
plot_borders(params, u_span, du_span, false)

%% Test #2
omega = -2;
L1 = 2 * pi / 3; L2 = pi / 3;
params = [omega L1 L2];

x_span_plus = [0, L1 + L2];
U_plus_set = f_scan(params, x_span_plus, u_span, du_span, M, N, intervals);

figure('Position', [100, 100, 350, 300])
axis([u_span du_span])
hold on

plot_scan(u_span, du_span, U_plus_set, false)
plot_borders(params, u_span, du_span, false)

%% Test #3
omega = -0.5;
L1 = pi / 3; L2 = 2 * pi / 3;
params = [omega L1 L2];

x_span_plus = [0, L1 + L2];
U_plus_set = f_scan(params, x_span_plus, u_span, du_span, M, N, intervals);

figure('Position', [100, 100, 350, 300])
axis([u_span du_span])
hold on

plot_scan(u_span, du_span, U_plus_set, false)
plot_borders(params, u_span, du_span, false)