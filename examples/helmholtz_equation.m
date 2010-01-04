% This is a script that tests DG solution of the 1D definite Helmholtz equation
%
%  k u - u_{xx} = f
%
% We do this using LDG: rewrite it as:
%
% k u - r_x = f
% r - u_x = 0
%
% The fluxes are then chosen as upwinding from opposite sides.
%
% The boundary conditions are periodic, though this is not necessary.

clear
close all

from labtools import spdiag
from labtools.linalg import block_spdiag
from labtools import typelatex

from dg.meshes import random_mesh_1d equidistant_mesh_1d local_operators
from dg.fluxes import right_value left_value flux_block

k = 1; % Should be positive for a definite problem
f = @(x) (k+1)*sin(x); % Forcing function
uexact = @(x) sin(x);
rexact = @(x) cos(x);

interval = [0, 2*pi];  % global x-interval 
K = 9;              % Number of cells
N = 4;              % (degree+1) of polynomial on each cell

% Step 1: create the mesh
mesh = equidistant_mesh_1d(interval, K, 'N', N);
%mesh = random_mesh_1d(interval, K, 'N', N);

% impose periodic boundary conditions:
mesh.face_to_face(1) = 2*K;  % "the exterior face on the left is the face on the right"
mesh.face_to_face(2*K) = 1;  % "the exterior face on the right is the face on the left"

% Step 2: Create mass, stiffness matrices on the mesh:
ops = local_operators(mesh.local_nodes);
jacobian = spdiag(1./mesh.cell_scale); % Scales for each element

% Step 3: Build operator 
%  This will be done using direct inversion of a matrix. To form the matrix,
%  we'll operate on unit vectors to construct each column.
%
%  Note that if we wanted to solve A*x = b using an indirect (iterative) method,
%  this step would be much easier than below ... but then step 4 is harder.

A = spalloc(2*K*N, 2*K*N, 2*N*2*K*N);  % About 2*N entries on each row

% Create some indices to make things easier later:
block1_start = 1;
block1_end = K*N;
block2_start = K*N+1;
block2_end = 2*K*N;
block1 = block1_start:block1_end;
block2 = block2_start:block2_end;

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

% rhs:
b = [f(mesh.nodes(:)); 0*mesh.nodes(:)];

% Step 4: solve system
x = A\b;

% plot solutions
figure; 
plot(mesh.nodes, reshape(x(block1), [N K]), '.'); hold on;
plot(mesh.nodes(:), uexact(mesh.nodes(:)), 'k--');
typelatex(title('Solution $u(x)$'));

figure; 
plot(mesh.nodes, reshape(x(block2), [N K]), '.'); hold on;
plot(mesh.nodes(:), rexact(mesh.nodes(:)), 'k--');
typelatex(title('Derivative $u''(x)$'));
