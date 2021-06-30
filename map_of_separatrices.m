%%
clc
clear

addpath('elliptic')

L1 = 2;
L2 = 1;
L = L1 + L2;

params = [-1 L1 L2];
u_span = [-15 +15];
u_sep = u_span(1) : 0.05 : u_span(2);

s_minus = @(u) (u / sqrt(2)) .* sqrt(2 + (u .^ 2));
s_plus = @(u) -s_minus(u);

du_sep_minus = s_minus(u_sep);
du_sep_plus  = s_plus(u_sep);

figure('Position', [100, 100, 350, 300])
hold on
plot(u_sep, du_sep_plus)
plot(u_sep, du_sep_minus)

n = length(u_sep);
u_map  = zeros(1, n);
du_map = zeros(1, n);

intervals = 2048;
x_span = [L1 L1 + L2];

for i = 1:n
	u0 = u_sep(i); du0 = du_sep_minus(i);
	[~, U] = f_solve(params, x_span, [u0 du0], intervals);
	u_map(i)  = U(end, 1);
	du_map(i) = U(end, 2);
end

plot(u_map, du_map)

a = 1;
k = 1 / sqrt(2);

n_plus = [0 1 2 3 4];
X_plus = icn(sqrt(a) / ((a^2 + 1) ^ (1/4)), k) + K(k) .* n_plus;
u_intersect_plus = (2 * a ^ (3/2)) / ((a^2 + 1) ^ (1/4)) .* X_plus / L2;

n_minus = [0 1 2 3 4];
X_minus = icn(sqrt(a) / ((a^2 + 1) ^ (1/4)), k) + K(k) .* n_minus;
u_intersect_minus = -(2 * a ^ (3/2)) / ((a^2 + 1) ^ (1/4)) .* X_minus / L2;

u_intersect = [u_intersect_minus u_intersect_plus];
du_intersect = s_plus(u_intersect);
plot(u_intersect, du_intersect, 'o')

axis([-12 +12 -90 +90])







