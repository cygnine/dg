function[mesh] = random_mesh_1d(interval,K, varargin)
% random_mesh_1d -- creates a random 1D finite element mesh
%
% mesh = random_mesh_1d(interval,K,{N=5,local_nodes=[]})
%
%   Creates a finite element mesh on the 1D interval 'interval' with K elements.
%   Each cell is assumed to have N degrees of freedom. The size of each element
%   is a random. All the necessary
%   connectivity information is stored in the struct mesh. 
%
%   If a set of local nodes on the standard interval [-1,1] is given
%   (local_nodes), this template is used to generate a global nodal set;
%   otherwise, the Legendre Gauss-Lobatto nodes with the given value of N are
%   used.

persistent glq strict_inputs mesh_from_cell_boundaries
if isempty(glq)
  from dg.meshes import mesh_from_cells_1d
  from labtools import strict_inputs
end

opt = strict_inputs({'local_nodes', 'N'}, {[], 5}, [], varargin{:});

cell_boundaries = (sort(rand([K-1 1]))-0.5)*diff(interval) + mean(interval);
cell_boundaries = [interval(1); cell_boundaries; interval(2)];
%cell_boundaries = linspace(interval(1), interval(2), K+1).';
cells = [cell_boundaries(1:K) cell_boundaries(2:end)];
if isempty(opt.local_nodes)
  mesh = mesh_from_cells_1d(cells, 'N', opt.N);
else
  mesh = mesh_from_cells_1d(cells, 'local_nodes', opt.local_nodes);
end
