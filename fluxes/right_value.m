function[fu] = right_value(u_minus, u_plus, normal_minus, normal_plus)
% right_value -- The right-hand side value as flux
%
% fu = right_value(u_minus, u_plus, normal_minus, normal_plus)
%
%     Returns the `right-hand' value as the flux. The normal vectors are needed
%     to determine which side is `left' and which is `right'. In this function,
%     `minus' and `plus' denote `interior' and `exterior', not `left' and
%     `right'.
%
%     This function only makes sense in one spatial dimension. The
%     implementation for this function is f(u) = {{u}} + [[u]].
%
%     Fluxes are to be used to approximate f(u) for the equation
%
%         u_t + f(u)_x = 0
%
%     Ex: f(u) = 1/2*u^2:
%
%         u_t + u u_x = 0

persistent u_jump u_mean
if isempty(u_jump)
  from dg.fluxes import u_jump u_mean
end

% The right-hand values have normal vectors + 1
fu = u_mean(u_minus, u_plus) + ...
     u_jump(j_minus, u_plus, normal_minus, normal_plus);
