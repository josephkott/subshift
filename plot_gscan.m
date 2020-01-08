function plot_gscan(params, u_span, du_span)

L1 = params(2); L2 = params(2);
L = L1 + L2;

intervals_scan = 512;
intervals_solve = 2048;

% M = 1024; N = 1024;
M = 128; N = 128; % debug mode

x_span_plus  = [0 +L];
x_span_minus = [0 -L];

U_plus_set = f_scan(params, x_span_plus, u_span, du_span, M, N, intervals_scan);
U_minus_set = f_scan(params, x_span_minus, u_span, du_span, M, N, intervals_scan);

result_image = create_image_from_scan(U_plus_set, U_minus_set);

% selector that specifies non-collapsing forward / backward points
selector = find(U_plus_set == 0 & U_minus_set == 0);
[du_indices, u_indices] = ind2sub([N M], selector);
pair = [du_indices u_indices];

% bind parameters
to_values  = @(u_index, du_index) map_indices_to_values(u_span, du_span, M, N, u_index, du_index);
to_indices = @(u, du) map_values_to_indices(u_span, du_span, M, N, u, du);

result_image_g1 = result_image;
result_image_g2 = result_image;

for i = 1:length(u_indices)
	du_index = pair(i, 1);
	u_index  = pair(i, 2);
	
	[u0, du0] = to_values(u_index, du_index);
	[~, U] = f_solve(params, x_span_minus, [u0 du0], intervals_solve);
	
	u_map  = U(end, 1);
	du_map = U(end, 2);
	
	[u_index_map, du_index_map] = to_indices(u_map, du_map);
	
	if ~isnan(u_index_map) && ~isnan(du_index_map)
		[g1_sign, g2_sign] = get_signs_of_g_functions(params, u_map, du_map);
		
		% Coloring
		if g1_sign == +1
			result_image_g1(du_index_map, u_index_map) = 0.75; % orange
		else
			result_image_g1(du_index_map, u_index_map) = 0.15; % ?
		end
		
		if g2_sign == +1
			result_image_g2(du_index_map, u_index_map) = 0.75; % orange
		else
			result_image_g2(du_index_map, u_index_map) = 0.15; % ?
		end
	end
end

f = figure('Position', [100, 100, 1500, 600]);
subplot(1, 2, 1)
plot_scan(u_span, du_span, result_image_g1, false, 'jet')
subplot(1, 2, 2)
plot_scan(u_span, du_span, result_image_g2, false, 'jet')

end

% ---------------------------------------------------------------------------------------------------- %

function [u, du] = map_indices_to_values(u_span, du_span, M, N, u_index, du_index)

if u_index <= 0
	disp 'ERROR: u_index invalid'
	return
end


if du_index <= 0
	disp 'ERROR: du_index invalid'
	return
end

if u_index > M || du_index > N
	u = NaN; du = NaN;
	return
end

u  =  u_span(1) + (( u_span(2) -  u_span(1)) / (M - 1)) * ( u_index - 1);
du = du_span(1) + ((du_span(2) - du_span(1)) / (N - 1)) * (du_index - 1);

end

function [u_index, du_index] = map_values_to_indices(u_span, du_span, M, N, u, du)

if u < u_span(1) || u > u_span(2)
	u_index = NaN; du_index = NaN;
	return
end

if du < du_span(1) || du > du_span(2)
	u_index = NaN; du_index = NaN;
	return
end

u_index  = round((M - 1) * ( u -  u_span(1)) / ( u_span(2) -  u_span(1))) + 1;
du_index = round((N - 1) * (du - du_span(1)) / (du_span(2) - du_span(1))) + 1;

end

function result_image = create_image_from_scan(U_plus_set, U_minus_set)

result_image = zeros(size(U_plus_set));

result_image(U_plus_set == 1 & U_minus_set == 1) = 0; % collapse!
result_image(U_plus_set == 0) = 0.3;
result_image(U_minus_set == 0) = 0.6;
result_image(U_plus_set == 0 & U_minus_set == 0) = 1; % intersection of non-collapsing forward / backward

end


function [g1_sign, g2_sign] = get_signs_of_g_functions(params, u0, du0)
% g_1(p) = (DT_p e_1, e_1) \dot (DT_p e_2, e_1)
% g_2(p) = (DT_p e_1, e_2) \dot (DT_p e_2, e_2)
% where p = (u_0, u_0')

L1 = params(2); L2 = params(2);
L = L1 + L2;

% More dense grid than during scan
intervals = 2048;

% Create I-st quadrant vector frame
delta = 1e-3;

p0 = [u0 du0];
p1 = [u0 + delta, du0];
p2 = [u0, du0 + delta];

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

g1_sign = sign(dot(e1_map, e1) * dot(e2_map, e1));
g2_sign = sign(dot(e1_map, e2) * dot(e2_map, e2));

end

% ---------------------------------------------------------------------------------------------------- %
