function u_end = map(x, u0, du0)

if (x == 0)
	u_end = [u0 du0];
	return
end

ode_problem = @(x, u) [u(2); u(1) - (u(1) .^ 3)];
% [~, U] = ode4(ode_problem, [0, x], [u0; du0], 2 ^ 12);
[~, U] = ode45(ode_problem, [0, x], [u0; du0]);

u_end = U(end, :);

end

