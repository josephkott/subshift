function u = dn(x, k)

n = length(k);
if n == 1
	[~, ~, u] = ellipj(x, k ^ 2);
else
	m = length(x);
	u = zeros(n, m);
	for i = 1:n
		[~, ~, u(i, :)] = ellipj(x, k(i) ^ 2);
	end
end

end