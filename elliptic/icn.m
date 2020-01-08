function x = icn(a, k)
% Solve cn(x, k) = a -> x = cn^{-1}(a, k)

f = @(t) 1 ./ sqrt(1 - (k ^ 2) * (sin(t) .^ 2));
x = integral(f, 0, acos(a));

end

