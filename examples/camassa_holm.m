% This is a script that tests DG evolution of the Camassa-Holm equation
%
%       q_t = p_x - f_x(u) - B_x(r),
%             r = u_x
%             p = (r*u)_x
%
% The boundary conditions are periodic.

clear
close all

from labtools import spdiag
from odesolve.coeffs import lserk4

from dg import plot_solution
from dg.meshes import random_mesh_1d equidistant_mesh_1d local_operators
from dg.drivers import camassa_holm_ldg as rhs_driver
from dg.examples import helmholtz_matrix

u0 = @(x) 0.25*exp(-abs(x));  % Initial data
kappa = 0;
uexact = @(x,t) 0.25*exp(-abs(x-0.25*t));

interval = [-25, 25];  % global x-interval 
K = 100;              % Number of cells
N = 5;              % (degree+1) of polynomial on each cell

% Step 1: create the mesh
mesh = equidistant_mesh_1d(interval, K, 'N', N);
%mesh = random_mesh_1d(interval, K, 'N', N);

% impose periodic boundary conditions:
mesh.face_to_face(1) = 2*K;  % "the exterior face on the left is the face on the right"
mesh.face_to_face(2*K) = 1;  % "the exterior face on the right is the face on the left"

% Must recompute normal vectors having changed face connections
mesh.normal_plus = mesh.face_normals(mesh.face_to_face);

% Step 2: Create mass, stiffness matrices on the mesh:
ops = local_operators(mesh.local_nodes);
jacobian = spdiag(1./mesh.cell_scale); % Scales for each element

% Now get Helmholtz solver:
A = helmholtz_matrix(1, mesh, ops, jacobian);
A = inv(A);
A = A(1:end/2, 1:end/2);

% Step 3: specify time evolution operator 
rk = lserk4();   % Gets R-K coefficients

t0 = 0;  % initial time
t = t0;  % current time
T = 3;   % terminal time

dt = 0.5*mesh.dx^2;  % This is a heuristic, but it turns out to be approximately
                            % correct. If you don't decrease dt as you increase N, you'll
                            % get an instability very quickly.

% Step 4: start solving
u = u0(mesh.nodes);  % Evaluate initial data
myplots = plot(mesh.nodes, u, '.'); hold on;
eplot = plot(mesh.nodes(:), uexact(mesh.nodes(:),t), 'k--');
axis_limits = [mesh.interval, min(u(:)), max(u(:))];

ku = zeros(size(u)); % Allocating storage for RK

pause

while t<T;
  if (t+dt)>T
    dt = (T-t);
  end

  % RK stages:
  for q = 1:rk.p
    stage_time = t + dt*rk.c(q);   % Unnecessary since rhs doesn't depend on t

    % Compute numerical RHS:
    rhs = rhs_driver(u, mesh, ops, jacobian, A, kappa);

    % Update RK data
    ku = rk.a(q)*ku + dt*rhs;
    u = u + rk.b(q)*ku;
  end

  t = t + dt;

  %%%%%%%%%%% Plotting
  set(eplot, 'ydata', uexact(mesh.nodes(:),t)); 
  plot_solution(myplots, u, axis_limits);
  %%%%%%%%%% End plotting
end
