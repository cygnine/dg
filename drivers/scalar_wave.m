function[rhs] = scalar_wave(u, mesh, ops, jacobian, pde_flux, numerical_flux)
% scalar_wave -- rhs driver for scalar wave equation
%
% rhs = scalar_wave(u, mesh, ops, jacobian, pde_flux, numerical_flux)
%
%     Computes the right-hand side of the DG approximation to the scalar wave
%     equation. The function pde_flux should take only one input, u. The
%     function numerical_flux should be a function of five variables: u_minus,
%     u_plus, normal_minus, normal_plus, and pde_flux.

u_minus = u(mesh.face_indices);       % "interior" evaluations
u_plus = u_minus(mesh.face_to_face);  % "exterior" evaluations

% Evaluate the flux:
flux = numerical_flux(u_minus, u_plus, mesh.normal_minus, mesh.normal_plus, pde_flux);
flux = reshape(flux, [2 mesh.K]);
u_minus = reshape(u_minus, [2, mesh.K]);

% Evaluate the "local" rhs:
rhs = -ops.strong_diffmat*(pde_flux(u)) + ops.liftmat*(pde_flux(u_minus) - flux); % "strong" form
%rhs = ops.weak_diffmat*(pde_flux(u)) - ops.liftmat*(flux); % "weak" form

rhs = rhs*jacobian;   % scale appropriately to "global" rhs
