% This is a script that tests convergence of the DG solution of the 1D definite
% Helmholtz equation
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

from labtools import spdiag typelatex
from speclab.orthopoly1d.jacobi.eval import eval_jacobi_poly
from speclab.orthopoly1d.jacobi.quad import gauss_lobatto_quadrature as glq

from dg.meshes import replicate_local_nodes as repnodes
from dg.meshes import equidistant_mesh_1d local_operators
from dg.fluxes import right_value left_value flux_block
from dg.examples import helmholtz_matrix

k = 1; % Should be positive for a definite problem
f = @(x) (k+1)*sin(x); % Forcing function
uexact = @(x) sin(x);
rexact = @(x) cos(x);

interval = [0, 2*pi];  % global x-interval 

Ks = [5:50];
Ns = [3:10];
Kmin = min(Ks);
Nmin = min(Ns);

L2_errors = zeros([length(Ks) length(Ns)]);
L1_errors = zeros([length(Ks) length(Ns)]);
Linf_errors = zeros([length(Ks) length(Ns)]);

% These are local interpolation nodes
[r,w] = glq(50);

for K = Ks;
  for N = Ns
    mesh = equidistant_mesh_1d(interval, K, 'N', N);

    % impose periodic boundary conditions:
    mesh.face_to_face(1) = 2*K;  % "the exterior face on the left is the face on the right"
    mesh.face_to_face(2*K) = 1;  % "the exterior face on the right is the face on the left"

    ops = local_operators(mesh.local_nodes);
    jacobian = spdiag(1./mesh.cell_scale); % Scales for each element

    % Create a local operator that interpolates to a more refined grid for computing
    % L_infinity errors.
    Vtemp = eval_jacobi_poly(r, 0:(N-1), 'alpha', 0, 'beta', 0);
    interpolate = Vtemp*ops.invV;
    R = repnodes(r, mesh.cells);

    A = helmholtz_matrix(k, mesh, ops, jacobian);

    % rhs:
    b = [f(mesh.nodes(:)); 0*mesh.nodes(:)];

    x = A\b;

    % Compute errors:
    u = reshape(x(1:(K*N)), [N K]);
    L2_errors(K-Kmin+1, N-Nmin+1) = sqrt(sum(sum((u-uexact(mesh.nodes)).^2.*mesh.weights)));
    L1_errors(K-Kmin+1, N-Nmin+1) = sum(sum(abs(u-uexact(mesh.nodes)).*mesh.weights));

    u = interpolate*u;
    Linf_errors(K-Kmin+1, N-Nmin+1) = max(abs(u(:) - uexact(R(:))));
  end
end

% Now plot errors:
% The following three figures are `classical' convergence plots: some error
% measurement vs number of degrees of freedom. For each fixed N (order of
% polynomial on each element) the plot on a log-log scale looks linear -- this
% means that there is convergence of some order p.  You see that as N is
% increased, the slope gets more vertical -- this means the order is increasing.
% This is the advantage of DG methods: both "h" adaptivity (travelling along one
% linear line) and "p" adaptivity (jumping from one line to the next). 
figure;
loglog(Ks, L2_errors);
temp = axis;
axis([min(Ks), max(Ks), 1e-16 temp(4)]);
typelatex(xlabel('Number of cells $K$'));
typelatex(ylabel('$L^2$ error'));
typelatex(title('$L^2$ error for increasing order of approximation $N$'));

figure;
loglog(Ks, L1_errors);
temp = axis;
axis([min(Ks), max(Ks), 1e-16 temp(4)]);
typelatex(xlabel('Number of cells $K$'));
typelatex(ylabel('$L^1$ error'));
typelatex(title('$L^1$ error for increasing order of approximation $N$'));

figure;
loglog(Ks, Linf_errors);
temp = axis;
axis([min(Ks), max(Ks), 1e-16 temp(4)]);
typelatex(xlabel('Number of cells $K$'));
typelatex(ylabel('$L^\infty$ error'));
typelatex(title('$L^\infty$ error for increasing order of approximation $N$'));
