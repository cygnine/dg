% This is a script that tests DG evolution of the 1D scalar wave equation
%
%  u_t + a * u_x = 0
%  u_t + (a*u)_x = u_t + (f(u))_x = 0
%
% where a is the constant wave speed. The boundary conditions are periodic.

clear
close all

from labtools import spdiag
from odesolve.coeffs import lserk4

from dg import plot_solution
from dg.meshes import random_mesh_1d equidistant_mesh_1d local_operators
from dg.fluxes import lax_friedrichs as lf

a = 2;  % wave speed
u0 = @(x) 1+ sin(2*pi*x);  % initial data
f = @(u) a*u;  % The flux function
uexact = @(x,t) 1 + sin(2*pi*(x-a*t));

interval = [0, 1];  % global x-interval 
K = 9;              % Number of cells
N = 6;              % (degree+1) of polynomial on each cell

% Note that all of the steps below are now (in principle) independent of the
% pde, save the choice of flux

% Step 1: create the mesh
%mesh = equidistant_mesh_1d(interval, K, 'N', N);
mesh = random_mesh_1d(interval, K, 'N', N);
% impose periodic boundary conditions:
mesh.face_to_face(1) = 2*K;  % "the exterior face on the left is the face on the right"
mesh.face_to_face(2*K) = 1;  % "the exterior face on the right is the face on the left"

% Extract the normal vectors:
normal_minus = mesh.face_normals;
normal_plus = mesh.face_normals(mesh.face_to_face);

% Step 2: Create mass, stiffness matrices on the mesh:

ops = local_operators(mesh.local_nodes);

% Note that since these are all *local* matrices, they operate correctly for functions
% defined on [-1,1]. In order to get the affine scale correct on each cell,
% we'll have to multiply by the jacobian, which we'll do on the fly in the code.

jacobian = spdiag(1./mesh.cell_scale); % Scales for each element

% Step 3: specify time evolution operator 
rk = lserk4();   % Gets R-K coefficients

t0 = 0;  % initial time
t = t0;  % current time
T = 3;   % terminal time

dt = 0.5*mesh.dx/(abs(a));  % This is a heuristic, but it turns out to be approximately
                            % correct. If you don't decrease dt as you increase N, you'll
                            % get an instability very quickly.

% Step 4: start solving
u = u0(mesh.nodes);  % Evaluate initial data
myplots = plot(mesh.nodes, u, '.'); hold on;
eplot = plot(mesh.nodes(:), uexact(mesh.nodes(:),t), 'k');
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

    % Gather facial evaluations:
    u_minus = u(mesh.face_indices);       % "interior" evaluations
    u_plus = u_minus(mesh.face_to_face);  % "exterior" evaluations

    % Evaluate the flux:
    flux = lf(u_minus, u_plus, normal_minus, normal_plus, f, abs(a));
    flux = reshape(flux, [2 K]);
    u_minus = reshape(u_minus, [2, K]);

    % Evaluate the "local" rhs:
    rhs = -ops.strong_diffmat*(f(u)) + ops.liftmat*(f(u_minus) - flux); % "strong" form
    %rhs = ops.weak_diffmat*(f(u)) - ops.liftmat*(flux); % "weak" form

    rhs = rhs*jacobian;   % scale appropriately to "global" rhs

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
