%%
f_phase = @(u, v) v.^2 - u.^2 - 0.5*u.^4;
u = -2:0.05:2;
v = -3:0.05:3;
[U, V] = meshgrid(u, v);

PH = f_phase(U, V);

figure('Position', [100, 100, 350, 300])
colormap(gray(1));
contour(U, V, PH, 'LineWidth', 1);

%%
f_phase = @(u, v) v.^2 - u.^2 + 0.5*u.^4;
u = -2:0.01:2;
v = -3:0.01:3;
[U, V] = meshgrid(u, v);

PH = f_phase(U, V);

figure('Position', [100, 100, 350, 300])
colormap(gray(1));
contour(U, V, PH, 'LineWidth', 1);