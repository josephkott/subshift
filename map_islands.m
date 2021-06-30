function P_islands = map_islands(params, u_span, du_span, M, N, U_plus_set, U_minus_set, map_distance)

intervals_solve = 2048;

% selector that specifies non-collapsing forward / backward points
selector = find(U_plus_set == 0 & U_minus_set == 0);
[du_indices, u_indices] = ind2sub([N M], selector);
pair = [du_indices u_indices];

P_islands = ones(N, M);

% bind parameters
to_values  = @(u_index, du_index) map_indices_to_values(u_span, du_span, M, N, u_index, du_index);
to_indices = @(u, du) map_values_to_indices(u_span, du_span, M, N, u, du);

for i = 1:length(u_indices)
	du_index = pair(i, 1);
	u_index  = pair(i, 2);
	
	[u0, du0] = to_values(u_index, du_index);
	[~, U] = f_solve(params, [0 map_distance], [u0 du0], intervals_solve);
	
	u_map  = U(end, 1);
	du_map = U(end, 2);
	
	[u_index_map, du_index_map] = to_indices(u_map, du_map);
	
	if ~isnan(u_index_map) && ~isnan(du_index_map)
		P_islands(du_index_map, u_index_map) = 0;
	end
end

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

% ---------------------------------------------------------------------------------------------------- %
