function plot_scan(u_span, du_span, U, varargin)

if length(varargin) == 1
	new_figure = varargin{1};
	if new_figure 
		figure('Position', [100, 100, 350, 300])
	end
end

image(u_span, du_span, U, 'CDataMapping','scaled')
set(gca, 'YDir', 'normal')
colormap('gray')

end

