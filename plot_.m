function plot_(params, u_span, du_span)

L1 = params(2); L2 = params(2);
L = L1 + L2;

intervals_scan = 512;

% M = 2048; N = 2048;
 M = 1024; N = 1024;
% M = 128; N = 128; % debug mode

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

result_image_signs = result_image;

for i = 1:length(u_indices)
	du_index = pair(i, 1);
	u_index  = pair(i, 2);
	
	[u0, du0] = to_values(u_index, du_index);
	same_signs = get_signs_of_linearization_matrix_values(params, u0, du0);
		
	% Coloring
	if same_signs == +1
		result_image_signs(du_index, u_index) = 0.75; % orange
	elseif same_signs == -1
		result_image_signs(du_index, u_index) = 0.5; % ?
	end
end

f = figure('Position', [100, 100, 750, 600]);
plot_scan(u_span, du_span, result_image_signs, false, 'jet')

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


function same_signs = get_signs_of_linearization_matrix_values(params, u0, du0)
% g_1(p) = (DT_p e_1, e_1) \dot (DT_p e_2, e_1)
% g_2(p) = (DT_p e_1, e_2) \dot (DT_p e_2, e_2)
% where p = (u_0, u_0')

L1 = params(2); L2 = params(3);
L = L1 + L2;

% More dense grid than during scan
intervals = 2048;

% Create I-st quadrant vector frame
delta = 1e-6;

p0 = [u0 du0];
p1 = [u0 + delta, du0];
p2 = [u0, du0 + delta];

[~, U] = f_solve(params, [0 -L], p0, intervals);
p0_map = U(end, :);

[~, U] = f_solve(params, [0 -L], p1, intervals);
p1_map = U(end, :);

[~, U] = f_solve(params, [0 -L], p2, intervals);
p2_map = U(end, :);

% Create vectors
e1 = p1 - p0;
e2 = p2 - p0;
e1_map = p1_map - p0_map;
e2_map = p2_map - p0_map;

A = [e1_map' e2_map'];

% if all(A(:) > 0)
if all(sign(A(:)) == [+1; -1; -1; +1;])
	same_signs = +1;
% elseif all(A(:) < 0)
elseif all(sign(A(:)) == [-1; +1; +1; -1;])
	same_signs = -1;
else
	same_signs = 0;
end

end

% ---------------------------------------------------------------------------------------------------- %
