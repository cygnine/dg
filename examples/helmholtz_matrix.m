function[A] = helmholtz_matrix(k, mesh, ops, jacobian)
% helmholtz_matrix -- Creates DG Helmholtz operator
%
% A = helmholtz_matrix(mesh, ops, jacobian)
%
%     Creates a sparse matrix that is the numerical representation for the operator
%
%         L[u] = -u_{xx} + k u 
%
%     when it is discretized using the discontinuous Galerkin method. LDG fluxes
%     are used and periodic boundary conditions are imposed. 

persistent block_spdiag flux_block left_value right_value
if isempty(block_spdiag)
  from labtools.linalg import block_spdiag
  from dg.fluxes import flux_block left_value right_value
end

K = mesh.K;
N = mesh.N;

A = spalloc(2*K*N, 2*K*N, 2*N*2*K*N);  % About 2*N entries on each row

% Create some indices to make things easier later:
block1_start = 1;
block1_end = K*N;
block2_start = K*N+1;
block2_end = 2*K*N;
block1 = block1_start:block1_end;
block2 = block2_start:block2_end;

% Rewrite it as:
%
% k u - r_x = f
% r - u_x = 0
%
%    [ u ]
%  A*[   ] = b
%    [ r ]

% Let's do r - u_x = 0 first:
flux_function = @(x) -x;

A(block2, block2) = speye(K*N);
A(block2, block1) = block_spdiag(-ops.strong_diffmat, 1./mesh.cell_scale);

% Do the flux via iteration:
block = flux_block(flux_function, left_value, mesh, ops, jacobian);
A(block2,block1) = A(block2,block1) + block;

% k u - r_x = f
% We only have to recompute the flux_block since we're using a different flux
% now:
block = flux_block(flux_function, right_value, mesh, ops, jacobian);

A(block1,block1) = k*A(block2,block2);
A(block1,block2) = block_spdiag(-ops.strong_diffmat, 1./mesh.cell_scale) + ...
                   block;
