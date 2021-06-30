%% Try MATLAB solvers
clc; clear
bvpfcn = @(x, y) [y(2) -y(1)];
bcfcn  = @(ya, yb) [ya(1) yb(2) - 1];
guess  = @(x) [sin(x) cos(x)];

xmesh = linspace(0, 1, 10);
solinit = bvpinit(xmesh, [0; 1]);
sol = bvp4c(bvpfcn, bcfcn, solinit);
plot(sol.x, sol.y, '-o')

%% Apply to a one-island periodic solution
clc; clear
L1 = 2; L2 = 1; L = L1 +  L2;
eta = @(x) sigma(x, L1, L2);

u0 = -1.341; du0 = 2.118;
bvpfcn = @(x, u) [u(2), u(1) + eta(x) * (u(1) .^ 3)];
bcfcn  = @(ua, ub) [ua(1) - u0, ub(2) - du0];

xmesh   = linspace(0, L, 1000);
solinit = bvpinit(xmesh, [u0; du0]);
sol = bvp4c(bvpfcn, bcfcn, solinit);

new_du0 = sol.y(2, 1);

bcfcn  = @(ua, ub) [ua(2) - du0, ub(1) - u0];

xmesh   = linspace(0, -L, 1000);
solinit = bvpinit(xmesh, [u0; du0]);
sol = bvp4c(bvpfcn, bcfcn, solinit);

new_u0 = sol.y(1, 1);
% 
% figure
% plot(sol.x, sol.y, '-o')

%% One iteration
u0 = sol.y(1, end);
du0 = sol.y(2, 1);

bcfcn  = @(ua, ub) [ua(1) - u0, ub(2) - du0];
solinit = bvpinit(xmesh, [u0; du0]);
sol = bvp4c(bvpfcn, bcfcn, solinit);

plot(sol.x, sol.y, '-o')

%% Iterations
n = 10;
ui = zeros(2, n);
ui(1, 1) = sol.y(1, end);
ui(2, 1) = sol.y(2, 1);

for i = 2:n
	u0 = ui(1, i - 1);
	du0 = ui(2, i - 1);
	
	bcfcn  = @(ua, ub) [ua(1) - u0, ub(2) - du0];
	solinit = bvpinit(xmesh, [u0 0]);
	sol = bvp4c(bvpfcn, bcfcn, solinit);
	
	ui(1, i) = sol.y(1, end);
	ui(2, i) = sol.y(2, 1);
end

