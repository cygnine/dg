function[ops] = local_operators(local_nodes, varargin)
% local_operators -- Creates local elementwise operators
%
% ops = local_operators(local_nodes, {basis='jacobi', alpha=0, beta=0})
%
%     Creates elementwise mass and differentiation operators for computations.

persistent strict_inputs eval_jacobi_poly
if isempty(strict_inputs)
  from labtools import strict_inputs
  from speclab.orthopoly1d.jacobi.eval import eval_jacobi_poly
end

inputs = {'basis', 'alpha', 'beta'};
defaults = {'jacobi', 0, 0};
opt = strict_inputs(inputs, defaults, [], varargin{:});

N = length(local_nodes);

switch opt.basis
case 'jacobi'
  V = eval_jacobi_poly(local_nodes, 0:(N-1), 'alpha', opt.alpha, 'beta', opt.beta);
  invV = eye(N)/V;

  % The mass matrix is then defined as 
  local_mass_inv = V*V.';

  % The differentiation matrix
  dps = eval_jacobi_poly(local_nodes, 0:(N-1), 'alpha', opt.alpha, 'beta', opt.beta, 'd', 1);

  local_diffmat = dps*invV;

  local_stiffness = inv(local_mass_inv)*local_diffmat;
  weak_diffmat = local_mass_inv*(local_stiffness.');

  strong_diffmat = local_diffmat;

  temp = zeros([N 2]);
  temp(1,1) = -1; temp(N,2) = 1;
  liftmat = local_mass_inv*temp;

  ops.V = V;
  ops.invV = invV;
  ops.mass = inv(local_mass_inv);
  ops.invmass = local_mass_inv;
  ops.weak_diffmat = weak_diffmat;
  ops.strong_diffmat = strong_diffmat;
  ops.stiffness = local_stiffness;
  ops.liftmat = liftmat;
end
