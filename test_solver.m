%% Nonlinear potential plot
x = -3 * pi : 0.01 : 3 * pi;
L1 = pi / 2; L2 = pi / 2;
u = sigma(x, L1, L2);

plot(x, u, 'LineWidth', 2)
grid on

%% Test MEX solver
omega = -1;
L1 = pi / 2; L2 = pi / 2;
du = @(x, u) [u(2) -omega * u(1) - sigma(x, L1, L2) .* (u(1) ^ 3)];

intervals = 16;

[X0, U0] = ode4(du, [0 pi], [0.1 0.1], intervals);
plot(X0, U0, 'color', 'black', 'LineWidth', 2)

[X, U] = f_solve([omega L1 L2], [0 pi], [0.1 0.1], intervals);

hold on
plot(X, U, 'color', 'red')

diff_u  = max(abs(U0(:, 1) - U(:, 1)));
diff_du = max(abs(U0(:, 2) - U(:, 2)));

title(sprintf('Diff: u: %g, du: %g', diff_u, diff_du));