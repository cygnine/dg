function[fu] = godunov(u_minus, u_plus, g, h)
% godunov -- Godunov flux
%
% fu = godunov(u_minus, u_plus, g, h)
%
%     Uses the Godunov flux on the values u_minus and u_plus to approximate
%     f(u).
%
%     The flux function is
%
%     fu = g(u) if u_minus <= u_plus
%        = h(u) if u_minus > u_plus
%
%     g(u) is the minimum of f(u) over (u_minus, u_plus), and h(u) is the maximum
%     of f(u) over (u_plus, u_minus). g and h are function handles of two
%     inputs, u_minus and u_plus.
%
%     This function is vectorized as long as g and h are appropriately
%     vectorized for multiple rows of u_minus, u_plus.
% 
%     Fluxes are to be used to approximate f(u) for the equation
%
%         u_t + f(u)_x = 0
%
%     Ex: f(u) = 1/2*u^2;
%
%     g = @(a,b) 1/2*min([a.^2, b.^2], [], 2);
%     h = @(a,b) 1/2*max([a.^2, b.^2], [], 2);

flags = (u_minus<=u_plus);
fu = zeros(size(u_minus));

fu(flags) = g(u_minus(flags), u_plus(flags));

flags = ~flags;
fu(flags) = h(u_minus(flags), u_plus(flags));
