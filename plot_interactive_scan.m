function plot_interactive_scan(params, u_span, du_span)

omega = params(1); L1 = params(2); L2 = params(2);
L = L1 + L2;

intervals = 512;
% M = 1024; N = 1024;
M = 128; N = 128; % debug mode

x_span_plus  = [0 +L];
x_span_minus = [0 -L];

U_plus_set = f_scan(params, x_span_plus, u_span, du_span, M, N, intervals);
U_minus_set = f_scan(params, x_span_minus, u_span, du_span, M, N, intervals);

U_set = 1.2 * U_plus_set + U_minus_set;

f = figure('Position', [100, 100, 1000, 900]);
plot_scan(u_span, du_span, U_set, false)
hold on

mouse_click_event = @(~, ~) mouse_click_handler(params);
set(f, 'WindowButtonDownFcn', mouse_click_event)

end

function mouse_click_handler(params)

axes_handle = gca;
point = get(axes_handle, 'CurrentPoint');
u_click  = point(1, 1);
du_click = point(1, 2);

% Create I-st quadrant vector frame
delta = 1e-3;
p0 = [u_click du_click];
p1 = [u_click + delta, du_click];
p2 = [u_click, du_click + delta];

L1 = params(2); L2 = params(3); L = L1 + L2;
intervals = 2048;

[~, U] = f_solve(params, [0 +L], p0, intervals);
p0_map = U(end, :);

[~, U] = f_solve(params, [0 +L], p1, intervals);
p1_map = U(end, :);

[~, U] = f_solve(params, [0 +L], p2, intervals);
p2_map = U(end, :);

% Create vectors
e1 = p1 - p0;
e2 = p2 - p0;
e1_map = p1_map - p0_map;
e2_map = p2_map - p0_map;

% Normalise vectors
fprintf('e1_map: %.2f, %.2f\n', [e1_map(1) e1_map(2)])
fprintf('e2_map: %.2f, %.2f\n', [e2_map(1) e2_map(2)])
fprintf('norm: %.2f, %.2f\n', [norm(e1_map) norm(e2_map)])

e1_norm = e1 / norm(e1);
e2_norm = e2 / norm(e2);
e1_map_norm = e1_map / norm(e1_map);
e2_map_norm = e2_map / norm(e2_map);

angle = acos(dot(e1_map_norm, e2_map_norm)) * 180 / pi;
fprintf('%3.2f\n', angle)

% Plot on a small inner figure
axes('Position', [0.7 0.7 0.2 0.2])
box on; hold on

quiver(0, 0, e1_norm(1), e1_norm(2), 'Color', 'black')
quiver(0, 0, e2_norm(1), e2_norm(2), 'Color', 'black')
quiver(0, 0, e1_map_norm(1), e1_map_norm(2), 'Color', 'red', 'LineWidth', 2)
quiver(0, 0, e2_map_norm(1), e2_map_norm(2), 'Color', 'red', 'LineWidth', 2)

hold off

end