%% asymptotic formala for the 1st interval (*)

L1 = pi / 2; L2 = pi / 2; L = L1 + L2;
params = [-1 pi/2 pi/2];
u_span = [-4 +4]; du_span = [-8 +8];
M = 256; N = 256;

U_plus_set = f_scan(params, [0, L1], u_span, du_span, M, N, 1024);
plot_scan(u_span, du_span, U_plus_set, true) % not create new figure

%%

% du0 = +2;
du0 = +6.5;

% u0_min = -1.5; u0_max = -0.5;
u0_min = -2.9; u0_max = -2.46;
ugrid = u0_min:0.01:u0_max;

% eps = (du0 ^ 2) - ugrid .^ 2  - 0.5 * (ugrid .^ 4);

eps = 0.01;


discrepancy = zeros(1, length(eps));
u = zeros(1, length(eps));
du = zeros(1, length(eps));

params_ = [-1.5 7 pi/2];
for i = 1:length(eps)
	u0 = ugrid(i);
	[~, U] = f_solve(params, [0 7], [u0 du0], 2048); % ODE solving
	
	u(i) = U(end, 1);
	du(i) = U(end, 2);
	
	discrepancy(i) = (u(i) / (sqrt(2) + sqrt((u(i) ^ 2) + 2))) - ...
		(sqrt(2) + sqrt((u0 ^ 2) + 2)) / (-u0) * (exp(L1) * eps(i) / 32);
end

figure
% plot(eps, discrepancy);
plot(eps, u)

%%
clc
clear

u0 = -3;
eps = 0.0001;
du0 = sqrt((u0 ^ 2) + 0.5 * (u0 ^ 4) + eps);

L = 9;
intervals = 10000;
[~, U] = f_solve([-1 L 0], [0 L], [u0 du0], intervals);

if length(U) < (intervals + 1)
	u = Inf;
	du = Inf;
else
	u = U(end, 1);
	du = U(end, 2);
end

z = ((sqrt(2) + sqrt((u0 ^ 2) + 2)) / (-u0)) * (exp(L) * eps / 32);
u_asy = 2 * sqrt(2) * z / (1 - (z ^ 2));

fprintf('u = %g, u_asy = %g\n', [u, u_asy]);

%%
clc
clear

u0 = -3;
eps = 0.6;
du0 = sqrt((u0 ^ 2) + 0.5 * (u0 ^ 4) + eps);

L = 3;

[~, U] = f_solve([-1 L L], [0 L], [u0 du0], 2048);
plot(U)

%%
clc
clear

u0 = 2;
du0 = sqrt((u0 ^ 2) + 0.5 * (u0 ^ 4));
H0 = (du0 ^ 2) - (u0 ^ 2) + 0.5 * (u0 ^ 4);

K0 = K(1 / sqrt(2));
% du = @(x, u) [u(2); u(1) - (u(1) .^ 3)];

u_end = -u0;
du_end = du0;

ode_problem = @(x, u) [u(2); u(1) - (u(1) .^ 3)];
discrepancy = @(x) sqrt(sum((map(x, u0, du0) - [u_end du_end]) .^ 2));

% initial
% fminsearch(discrepancy, 0.5)
x0 = 2.220556640624998;

u0_grid = u0:0.1:6;
du0_grid = sqrt((u0_grid .^ 2) + 0.5 * (u0_grid .^ 4));
H0_grid = (du0 ^ 2) - (u0_grid .^ 2) + 0.5 * (u0_grid .^ 4);
x0_grid = zeros(1, length(u0_grid));
x0_grid(1) = x0;

for i = 2:length(u0_grid)
	u_end = -u0_grid(i);
	du_end = du0_grid(i);
	
	discrepancy = @(x) sqrt(sum((map(x, u0_grid(i), du0_grid(i)) - [u_end du_end]) .^ 2));
	x0_grid(i) = fminsearch(discrepancy, x0_grid(i - 1));
end

x0_asy = (2 * (2 ^ (3 / 4)) * K0 - 2.13006) ./ (H0_grid .^ (1 / 4));

%% Full

clc
clear

u0 = 1.5;
du0 = sqrt((u0 ^ 2) + 0.5 * (u0 ^ 4));
H0 = (du0 ^ 2) - (u0 ^ 2) + 0.5 * (u0 ^ 4);

K0 = K(1 / sqrt(2));
% du = @(x, u) [u(2); u(1) - (u(1) .^ 3)];

u_end = u0;
du_end = du0;

ode_problem = @(x, u) [u(2); u(1) - (u(1) .^ 3)];
discrepancy = @(x) sqrt(sum((map(x, u0, du0) - [u_end du_end]) .^ 2));

% initial
% fminsearch(discrepancy, 4)
x0 = 4.396972656250000;

u0_grid = u0:0.1:6;
H0_grid = (du0 ^ 2) - (u0_grid .^ 2) + 0.5 * (u0_grid .^ 4);
x0_grid = zeros(1, length(u0_grid));
x0_grid(1) = x0;

for i = 2:length(u0_grid)
	u_end = u0_grid(i);
	du_end = du0;
	
	discrepancy = @(x) sqrt(sum((map(x, u0_grid(i), du0) - [u_end du_end]) .^ 2));
	x0_grid(i) = fminsearch(discrepancy, x0_grid(i - 1), optimset('TolX', 1e-9));
end

x0_asy = (4 * K0) ./ ((2 * H0_grid) .^ (1 / 4));

%% Check the constant
int = @(x) sqrt(2) ./ sqrt(2 - (x .^ 4));
xspan = [-1 +1];

%% Another direction
clc
clear

u0 = -1.5;
du0 = sqrt((u0 ^ 2) + 0.5 * (u0 ^ 4));
H0 = (du0 ^ 2) - (u0 ^ 2) + 0.5 * (u0 ^ 4);

K0 = K(1 / sqrt(2));
% du = @(x, u) [u(2); u(1) - (u(1) .^ 3)];

u_end = -u0;
du_end = du0;

ode_problem = @(x, u) [u(2); u(1) - (u(1) .^ 3)];
discrepancy = @(x) sqrt(sum((map(x, u0, du0) - [u_end du_end]) .^ 2));

% initial
% fminsearch(discrepancy, 0.5)
x0 = 1.303027343750002;

u0_grid = u0:(-0.1):-6;
du0_grid = sqrt((u0_grid .^ 2) + 0.5 * (u0_grid .^ 4));
H0_grid = (du0_grid .^ 2) - (u0_grid .^ 2) + 0.5 * (u0_grid .^ 4);
x0_grid = zeros(1, length(u0_grid));
x0_grid(1) = x0;

for i = 2:length(u0_grid)
	u_end = -u0_grid(i);
	du_end = du0_grid(i);
	
	discrepancy = @(x) sqrt(sum((map(x, u0_grid(i), du0_grid(i)) - [u_end du_end]) .^ 2));
	x0_grid(i) = fminsearch(discrepancy, x0_grid(i - 1));
end

%for i = 1:length(u0_grid)
%	int = @(x) 1 ./ sqrt((x .^ 2) - 0.5 * (x .^ 4) + H0_grid(i));
%	
%	x0_asy(i) = integral(int, -((H0_grid(i))^(1/4)), +(H0_grid(i)^(1/4)));
%end

x0_asy = 2.13006 ./ ((H0_grid) .^ (1 / 4));

plot(H0_grid, x0_grid, H0_grid, x0_asy)