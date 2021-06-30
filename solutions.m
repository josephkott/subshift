%%
clc
clear

L1 = 2;
L2 = 1;
L = L1 + L2;
params = [-1 L1 L2];

xspan = [-5 * L,  -L2 / 2];
% xspan_right = [+2 * L + L1, 0];

c = 0:0.01:80;
left = zeros(length(c), 2);
% right = zeros(length(c), 2);

for i = 1:length(c)
	u0 = c(i) * exp(xspan(1));
	du0 = c(i) * exp(xspan(1));
	[X, U] = f_solve(params, xspan, [u0 du0], 8192);
	left(i, :) = U(end, :);
end

plot(c, left(:, 1))

% c = 0:0.01:7;
% for i = 1:length(c)
% 	u0 = c(i) * exp(-xspan_right(1));
% 	du0 = c(i) * exp(-xspan_right(1));
% 	[X, U] = f_solve(params, xspan_right, [u0 du0], 1024);
% 	right(i, :) = U(end, :);
% end
% 
% plot(left(:, 1), left(:, 2), right(:, 1), right(:, 2))

%%

f = @(c) get_ux_end(params, xspan, c);

eps = 1e-9;
c_left = 6.7;
c_right = 7.2;
c_root = dichotomy(f, c_left, c_right, eps);

%%

f = @(c) get_u_end(params, xspan, c);

eps = 1e-9;
c_left = 66.89;
c_right = 67.98;
c_root = dichotomy(f, c_left, c_right, eps);

%% solution

u0 = c_root * exp(xspan(1));
du0 = c_root * exp(xspan(1));
[X, U] = f_solve(params, [xspan(1) xspan(2)], [u0 du0], 8192);

figure('Position', [100, 100, 350, 300])
plot([X, -X(end:-1:1) - L2], [U(:, 1) U(end:-1:1, 1)], 'Color', 'black')
axis([-2*L +2*L 0 2.1])

xticks([-3*L -2*L -L 0 L 2*L 3*L])
xticklabels({'-3L','-2L','-L','0','+L','+2L','+3L'})

%%


u0 = c_root * exp(xspan(1));
du0 = c_root * exp(xspan(1));
[X, U] = f_solve(params, [xspan(1) xspan(2)], [u0 du0], 8192);

figure('Position', [100, 100, 350, 300])
plot([X, -X(end:-1:1) - L2], [U(:, 1) -U(end:-1:1, 1)], 'Color', 'black')
axis([-2*L +2*L -6 6])

xticks([-3*L -2*L -L 0 L 2*L 3*L])
xticklabels({'-3L','-2L','-L','0','+L','+2L','+3L'})
%% ticks example
x = linspace(-10,10,200);
y = cos(x);
plot(x,y)
xticks([-3*pi -2*pi -pi 0 pi 2*pi 3*pi])
xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
yticks([-1 -0.8 -0.2 0 0.2 0.8 1])
