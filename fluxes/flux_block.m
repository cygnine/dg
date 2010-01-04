function[block] = flux_block(pde_flux, flux_function, mesh, ops, jacobian, varargin)
% flux_block -- creates the (linear) flux operator matrix
%
% block = flux_block(pde_flux, flux_function, mesh, ops, jacobian, {type='strong'})
%
%     Returns a matrix that acts as the global flux operator for the strong/weak
%     nodal DG methods. This function is mainly used for constructing elliptic
%     operators. As of now, this function only works in 1D.
%
%     The function handle 'flux_function' assumes that only five inputs are needed:
%     u_minus, u_plus, normal_minus, normal_plus, and flux_function. 

K = mesh.K;
N = mesh.N;

% Extract the normal vectors:
normal_minus = mesh.face_normals;
normal_plus = mesh.face_normals(mesh.face_to_face);

block = spalloc(K*N,K*N,2*K*N);

% Do the flux via iteration:
for q = 1:K;
  u_minus = zeros([length(mesh.face_indices) 1]);

  % The left face:
  index = 1 + (q-1)*N;
  u_minus(1+2*(q-1)) = 1;
  u_plus = u_minus(mesh.face_to_face);
  flux = flux_function(u_minus, u_plus, normal_minus, normal_plus, pde_flux);
  flux = reshape(flux, [2 K]); u_minus = reshape(u_minus, [2 K]);
  temp = ops.liftmat*(flux - pde_flux(u_minus))*jacobian;
  block(:,index) = temp(:);

  % The right face:
  index = q*N;
  u_minus = zeros([length(mesh.face_indices) 1]);
  u_minus(2*q) = 1;
  u_plus = u_minus(mesh.face_to_face);
  flux = flux_function(u_minus, u_plus, normal_minus, normal_plus, pde_flux);
  flux = reshape(flux, [2 K]); u_minus = reshape(u_minus, [2 K]);
  temp = ops.liftmat*(flux - pde_flux(u_minus))*jacobian;
  block(:,index) = temp(:);
end
