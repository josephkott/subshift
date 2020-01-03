function plot_separatrix_f_map(params, u_span)

% TODO: use omega
omega = params(1); L1 = params(2); L2 = params(3);

u_step = 0.01;
u_sep = u_span(1) : u_step : u_span(2);

s_plus = @(u) (u / sqrt(2)) .* sqrt(2 + (u .^ 2));
s_minus = @(u) -s_plus(u);
du_sep = s_plus(u_sep);

n = length(u_sep);
u_map  = zeros(1, n);
du_map = zeros(1, n);

intervals = 2048;
x_span = [L1 L1 + L2];

for i = 1 : n
	% TODO: replace with f_map
	u0 = u_sep(i); du0 = du_sep(i);
	[~, U] = f_solve(params, x_span, [u0 du0], intervals);
	u_map(i)  = U(end, 1);
	du_map(i) = U(end, 2);
end

% Compute predicted separatrix map shape
k = 1 / sqrt(2);
delta = icn(-1 / nthroot(2, 4), k);
u_predicted = -nthroot(2, 4) * u_sep .* cn(delta + nthroot(2, 4) * L2 * u_sep, k);
du_predicted = sqrt(2) * (u_sep .^ 2) .* dn(delta + nthroot(2, 4) * L2 * u_sep, k) ...
	.* sn(delta + nthroot(2, 4) * L2 * u_sep, k);

% Compute predicted ponts of intersections s^{-}(u) with T_0 [ s^{+}(u) ]
m = [1 2 3 4 5];
u_intersect = (-2 * icn(-1 / nthroot(2, 4), k) + 4 * K(k) * m) / (nthroot(2, 4) * L2);
du_intersect = s_minus(u_intersect);

figure('Position', [100, 100, 350, 300])
hold on

plot(u_sep, s_minus(u_sep), '-', 'Color', 'black', 'LineWidth', 2)
plot(u_sep, du_sep, '--', 'Color', 'black')
plot(u_map, du_map, 'Color', 'blue')
plot(u_predicted, du_predicted, 'Color', 'red')
plot(u_intersect, du_intersect, '*', 'Color', 'red')

% TODO: same things for u_sep < 0 -> complete map of the separatrix on [L1, L1 + L2];

end

