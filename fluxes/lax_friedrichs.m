function[fu] = lf(u_minus, u_plus, normal_minus, normal_plus, f,alpha)
% lax_friedrichs -- The Lax-Friedrichs flux
%
% fu = lf(u_minus, u_plus, normal_minus, normal_plus, f, alpha)
%
%     Uses the Lax-Friedrichs function to evaluate \hat{f}(u) given u_minus,
%     u_plus, and a function handle f.
%
%     The flux function is
%
%     fu = {{f(u)}} + alpha/2 * [[u]],
%
%     where {{.}} is the average, and [[u]] is the jump. (The jump function
%     requires knowledge of the normal vectors on the grid faces.) The constant alpha
%     gives the maximum of |f'(u)| over the relevant values of u (I.e. the max
%     Jacobian).
%
%     This function is vectorized for multiple rows of u_minus and u_plus as
%     long as f is vectorized accordingly.
%
%     Fluxes are to be used to approximate f(u) for the equation
%
%         u_t + f(u)_x = 0
%
%     Ex: f(u) = 1/2*u^2;
%
%     f = @(x) 1/2*x.^2;
%     alpha = 1 (if, say u \in [-1,1]).

persistent u_jump u_mean
if isempty(u_jump)
  from dg.fluxes import u_jump u_mean
end

fum = f(u_minus);
fup = f(u_plus);

fu = u_mean(fum,fup) + ...
     alpha/2*u_jump(u_minus, u_plus , normal_minus, normal_plus);
