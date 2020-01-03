function s = sigma(x, L1, L2)

% Total period
L = L1 + L2;

n = floor(x / L);
x0 = x - n * L;

s = ones(1, length(x));
s(x0 < L1) = -1;

end

