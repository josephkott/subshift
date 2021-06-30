function m = dichotomy (func, a, b, eps)

fa = func(a);
fb = func(b);
while (((b-a)/2) >= eps)
	c = (a + b) / 2;
	fc = func(c);
	if sign(fc) ~= sign(fa)
		b = c;
		fb = fc;
	else
		a = c;
		fa = fc;
	end
end

m = (a + b) / 2;

end