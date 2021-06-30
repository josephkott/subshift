function  u_end = get_u_end(param, xspan, c)

u0 = [c*exp(xspan(1)); c*exp(xspan(1))];
[~,U] = f_solve(param, xspan, u0, 8192);
u_end = U(end, 1);

end

