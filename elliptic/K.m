function result = K(k, varargin)
% Vectorized version of ellipk

if length(varargin) == 1
	tol = varargin{1};
else
	tol = 2.22e-16; % default from help
end

n = length(k);
result = zeros(1, n);

for i = 1:n
	result(i) = ellipk(k(i), tol);
end

end

