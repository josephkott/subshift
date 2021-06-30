function  ux_end = get_ux_end(params, xspan, c)

u0 = [c*exp(xspan(1)); c*exp(xspan(1))];
[~,U] = f_solve(params, xspan, u0, 8192);
ux_end = U(end, 2);

end

