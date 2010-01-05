function[rhs] = camassa_holm_ldg(u, mesh, ops, jacobian, A, kappa)
% camassa_holm_ldg -- RHS driver for the Camassa-Holm equation
%
% rhs = camassa_holm_ldg(u, mesh, ops, jacobian, A, kappa, pde_flux)
%
%     The right-hand side driver for the equation:
%
%       q_t = p_x - f_x(u) - B_x(r),
%             r = u_x
%             p = (r*u)_x
%
%     The function f is defined as 
%
%         f(u) = 2*kappa*u + 3/2*u^2
%
%     The numerical fluxes are chosen as given in [1]. The matrix A connects the
%     independent variable u with the derived variable q. u and q are related by
%
%       -u_xx + u = q
%
%     The matrix A should represent the inverse of the LDG discretization of the
%     left-hand side. I.e., u = A*q.
% 
%   [1]: "A local discontinuous Galerkin method for the Camassa-Holm equation",
%        Y. Xu & C.-W. Shu

persistent right left lf
if isempty(right)
  from dg.fluxes import right_value as right
  from dg.fluxes import left_value as left
  from dg.fluxes import lax_friedrichs as lf
end

id = @(x) x;

u_minus = u(mesh.face_indices);
u_plus = u_minus(mesh.face_to_face);
flux = right(u_minus, u_plus, mesh.normal_minus, mesh.normal_plus, id);
flux = reshape(flux, [2, mesh.K]);
u_minus = reshape(u_minus, [2, mesh.K]);

% Compute r: 
% r - u_x = 0
% \hat{u} = u_{from the right}
r = ops.strong_diffmat*u + ops.liftmat*(flux - u_minus);
r = r*jacobian;

% Compute p via collocation:
% p - (b(r) * u)_x = 0, 
%
% where b(r) = r = 1/2 * d/dr (r^2) = B'(r)
%
% \hat{u} = u_{from the right}
% \hat{b} = \frac{B(r+) - B(r-)}{r+ - r-}
B = @(r) 1/2*r.^2;

r_minus = r(mesh.face_indices);
r_plus = r_minus(mesh.face_to_face);
r_right = right(r_minus, r_plus, mesh.normal_minus, mesh.normal_plus, id);
r_left = left(r_minus, r_plus, mesh.normal_minus, mesh.normal_plus, id);

b_minus = reshape(r_minus, [2, mesh.K]);
b_flux = (B(r_right) - B(r_left))./(r_right - r_left);
b_flux = reshape(b_flux, [2, mesh.K]);

% since the u-fluxes are the same, we don't have to recompute that
p = ops.strong_diffmat*(r.*u) + ops.liftmat*(b_flux.*flux - b_minus.*u_minus);
p = p*jacobian;

% Now evolution of q := (u + u_{xx})
% q_t = p_x - B_x(r) - f_x(u),
%
%   where B(r) = 1/2*r^2 and f(x) = 2*kappa*u + 3/2*u^2. Again, we use
%   collocation for the whole thing.

f = @(u) 2*kappa*u + 3/2*u.^2;
lf_flux = @(um, up, nm, np, f) lf(um, up, nm, np, f, 2*kappa+3*max(abs(um)));

fu = f(u);
Br = B(r);

u_minus = u_minus(:);
f_minus = f(u_minus);
B_minus = B(r_minus);
p_minus = p(mesh.face_indices);
p_plus = p_minus(mesh.face_to_face);

f_flux = lf_flux(u_minus, u_plus, mesh.normal_minus, mesh.normal_plus, f);
B_flux = B(r_left);
p_flux = left(p_minus, p_plus, mesh.normal_minus, mesh.normal_plus, id);

f_diff = reshape(f_flux - f_minus, [2, mesh.K]);
p_diff = reshape(p_flux - p_minus, [2, mesh.K]);
B_diff = reshape(B_flux - B_minus, [2, mesh.K]);

rhs = ops.strong_diffmat*(p - Br - fu) + ops.liftmat*(...
            p_diff - B_diff - f_diff);
rhs = rhs*jacobian;

% Finally, transfer q to u via the matrix inverse A:
rhs = reshape(A*rhs(:), size(rhs));
