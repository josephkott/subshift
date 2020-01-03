function [ T, Y ] = ode4( ode, tspan, y0, N )

n = length(y0);
h = (tspan(2) - tspan(1)) / N;
T = tspan(1):h:tspan(2);

Y = zeros(N+1,n);
k = zeros(4,n);

Y(1,:) = y0;
for i = 1:N
	k(1,:) = h*ode(T(i), Y(i,:));
	k(2,:) = h*ode(T(i) + h/2, Y(i,:) + 1/2 * k(1,:));
	k(3,:) = h*ode(T(i) + h/2, Y(i,:) + 1/2 * k(2,:));
	k(4,:) = h*ode(T(i) + h, Y(i,:) + k(3,:));
	Y(i+1,:) = Y(i,:) + (1/6)*(k(1,:) + 2*k(2,:) + 2*k(3,:) + k(4,:));
end

end

