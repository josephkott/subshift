function plot_borders(params, u_span, du_span, varargin)

omega = params(1); L1 = params(2);

% Now I don't have suitable convetion about sign of \omega parameter in equation.
% Here I use the following form u_{xx} - \omega u - u^3 = 0, \omega > 0.
omega_negation = -omega;

[u_I,   du_I  ] = area_I  (omega_negation, L1); 
[u_II,  du_II ] = area_II (omega_negation, L1);
[u_III, du_III] = area_III(omega_negation, L1);

u =  [ u_III   u_II   u_I];
du = [du_III, du_II, du_I];

% Filtering
sel = (u < u_span(1)) | (u > u_span(2));
u  =  u(~sel);
du = du(~sel);

sel = (du < du_span(1)) | (du > du_span(2));
u  =  u(~sel);
du = du(~sel);

u_symmetric = -u(end : -1 : 1);
du_symmetric = -du(end : -1 : 1);

if length(varargin) == 1
	new_figure = varargin{1};
	if new_figure 
		figure('Position', [100, 100, 350, 300])
	end
end

plot(u, du, u_symmetric, du_symmetric, 'LineWidth', 2);

end


function [u, du] = area_I(omega, L1)
% L1: critical value, first period on nonlinear potential

x = sqrt(omega) * L1;

theta = @(k) -x / sqrt(2 * (k ^ 2) - 1);

% Cirtical value of k
options = optimset('Display','off');
k_critical = fsolve(@(k) theta(k) + 2 * K(k), (1 / sqrt(2)) + 0.1, options);

k_step = -0.0001;
k = 1 + k_step : k_step : k_critical - k_step;

% Dimensionless form
u_dless = zeros(1, length(k));
du_dless = zeros(1, length(k));

for i = 1 : length(k)
	u_dless(i) = - (sqrt(2) / sqrt(2 * (k(i) ^ 2) - 1)) ...
		* dn(theta(k(i)), k(i)) ...
		/ sn(theta(k(i)), k(i));
	du_dless(i) = (sqrt(2) / (2 * (k(i) ^ 2) - 1)) ...
		* cn(theta(k(i)), k(i)) ...
		/ (sn(theta(k(i)), k(i)) ^ 2);
end

u = sqrt(omega) * u_dless;
du = omega * du_dless;

end


function [u, du] = area_II(omega, L1)
% L1: critical value, first period on nonlinear potential

x = sqrt(omega) * L1;

theta = @(k) (sqrt(2 * (k ^ 2) + 2) * K(k) - x) / sqrt(2 * (k ^ 2) + 2);

k_critical = 1;
k_step = 0.0001;
k = 0 + k_step : k_step : k_critical - k_step;

% Dimensionless form
u_dless = zeros(1, length(k));
du_dless = zeros(1, length(k));

for i = 1 : length(k)
	u_dless(i) = ((1 - (k(i) ^ 2)) / sqrt((k(i) ^ 2) + 1)) ...
		* sn(theta(k(i)), k(i)) ...
		/ (cn(theta(k(i)), k(i)) * dn(theta(k(i)), k(i)));
	du_dless(i) = (sqrt(2) / (2 * (k(i) ^ 2) + 2)) ...
		* (((k(i) ^ 2) * (sn(theta(k(i)), k(i)) ^ 4) - 1) * ((k(i) ^ 2) - 1)) ...
		/ ((cn(theta(k(i)), k(i)) ^ 2) * (dn(theta(k(i)), k(i)) ^ 2));
end

u = sqrt(omega) * u_dless;
du = omega * du_dless;

end


function [u, du] = area_III(omega, L1)
% L1: critical value, first period on nonlinear potential

% Dimensionless form
x = sqrt(omega) * L1;

theta = @(k) K(k) - (x / sqrt(2 - 4 * (k ^ 2)));

% Cirtical value of k
options = optimset('Display','off');
k_critical = fsolve(@(k) theta(k) + K(k), 1 / (1.1 * sqrt(2)), options);

k_step = -0.0001;
k = k_critical + k_step : k_step : 0 - k_step;

% Dimensionless form
u_dless = zeros(1, length(k));
du_dless = zeros(1, length(k));

for i = 1 : length(k)
	u_dless(i) = 1 / sqrt(1 - (2 * (k(i) ^ 2))) ...
		* (sn(theta(k(i)), k(i)) * dn(theta(k(i)), k(i))) ...
		/ cn(theta(k(i)), k(i));
	du_dless(i) = (sqrt(2) / (2 - 4 * (k(i) ^ 2))) ...
		* ( (k(i) ^ 2) * (sn(theta(k(i)), k(i)) ^ 4) ...
			- 2 * (k(i) ^ 2) * (sn(theta(k(i)), k(i)) ^ 2) ...
			+ 1 ) ...
		/ (cn(theta(k(i)), k(i)) ^ 2);
end

u = sqrt(omega) * u_dless;
du = omega * du_dless;

end

