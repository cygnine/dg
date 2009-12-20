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

from speclab.orthopoly1d.jacobi.eval import eval_jacobi_poly

from dg import plot_solution
from dg.meshes import equidistant_mesh_1d
from dg.fluxes import lax_friedrichs as lf

a = 2;  % wave speed
u0 = @(x) 1+ sin(2*pi*x);  % initial data
f = @(u) a*u;  % The flux function
uexact = @(x,t) 1 + sin(2*pi*(x-a*t));

interval = [0, 1];  % global x-interval 
K = 9;              % Number of cells
N = 4;              % (degree+1) of polynomial on each cell

% Note that all of the steps below are now (in principle) independent of the
% pde, save the choice of flux

% Step 1: create the mesh
mesh = equidistant_mesh_1d(interval, K, 'N', N);
% impose periodic boundary conditions:
mesh.face_to_face(1) = 2*K;  % "the exterior face on the left is the face on the right"
mesh.face_to_face(2*K) = 1;  % "the exterior face on the right is the face on the left"

% Extract the normal vectors:
normal_minus = mesh.face_normals;
normal_plus = mesh.face_normals(mesh.face_to_face);

% Step 2: Create mass, stiffness matrices on the mesh:

% Form the local Vandermonde matrix:
V = eval_jacobi_poly(mesh.local_nodes, 0:(N-1), 'alpha', 0, 'beta', 0);

% The matrix V takes point evaluations of a function and transforms them into
% the expansions coefficients of a Legendre polynomial expansion. To get the
% inverse transformation:
invV = eye(N)/V;  % same as inv(V)

% The mass matrix is then defined as 
local_mass_inv = V*V.';

% The differentiation matrix:
% This matrix takes point evaluations to expansion coefficients, differentiates
% the expansion, then interpolates back to the point locations.
dps = eval_jacobi_poly(mesh.local_nodes, 0:(N-1), 'alpha', 0, 'beta', 0, 'd', 1);
local_diffmat = dps*invV;

local_stiffness = inv(local_mass_inv)*local_diffmat;
weak_diffmat = local_mass_inv*(local_stiffness.');

% Note that since this is a *local* matrix, it operates correctly for functions
% defined on [-1,1]. In order to get the affine scale correct on each cell,
% we'll have to multiply by the jacobian, which we'll do on the fly in the code.

jacobian = spdiag(1./mesh.cell_scale); % Scales for each element

% Last thing: must compute a matrix to "lift" boundary fluxes to interior of
% elements. It acts as the boundary integral in ibp. It's a local matrix: it
% operates on one element at a time.
temp = zeros([N 2]); 
temp(1,1) = -1; temp(N,2) = 1;  
liftmat = local_mass_inv*temp;

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
    rhs = -local_diffmat*(f(u)) + liftmat*(f(u_minus) - flux); % "strong" form
    %rhs = weak_diffmat*(f(u)) - liftmat*(flux); % "weak" form

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
