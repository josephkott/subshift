function plot_scan(u_span, du_span, U, varargin)

if length(varargin) >= 1
	new_figure = varargin{1};
	if new_figure 
		figure('Position', [100, 100, 350, 300])
	end
end

if length(varargin) >= 2
	specified_colormap = varargin{2};
else
	specified_colormap = 'jet';
end

image(u_span, du_span, U, 'CDataMapping','scaled')
set(gca, 'YDir', 'normal')
colormap(specified_colormap)

end

