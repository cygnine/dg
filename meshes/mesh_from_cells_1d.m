function[mesh] = mesh_from_cells_1d(cells, varargin)
% mesh_from_cells_1d -- The main 1D mesh constructor
%
% mesh = mesh_from_cells_1d(cells, {local_ndoes=[], N=5});
%
%   Creates a 1D finite element mesh with N degrees of freedom on each cell,
%   where the cells are defined by cells.  cells should be a K x 2 matrix,
%   meaning K cells with each row delineating the interval of the cell. All the
%   necessary connectivity information is stored in the struct mesh. 
%
%   If a set of local nodes on the standard interval [-1,1] is given
%   (local_nodes), this template is used to generate a global nodal set;
%   otherwise, the Legendre Gauss-Lobatto nodes are used.

persistent repnodes glq strict_inputs spdiag
if isempty(repnodes)
  from dg.meshes import replicate_local_nodes as repnodes
  from speclab.orthopoly1d.jacobi.quad import gauss_lobatto_quadrature as glq
  from labtools import strict_inputs spdiag
end

opt = strict_inputs({'local_nodes','local_weights', 'N'}, {[],[],5}, [], varargin{:});

if isempty(opt.local_nodes)
  [opt.local_nodes, opt.local_weights] = glq(opt.N, 'alpha', 0, 'beta', 0);
else
  opt.N = length(opt.local_nodes);
end

mesh.N = opt.N;
mesh.K = size(cells, 1);
mesh.cells = cells;
mesh.interval = [min(min(cells)), max(max(cells))];
mesh.cell_boundaries = sort([cells(:,1); mesh.interval(2)]);

mesh.cell_scale = diff(mesh.cells,1,2)/2;
mesh.cell_shift = mean(mesh.cells, 2);

mesh.nodes = repnodes(opt.local_nodes, mesh.cells);
mesh.weights = repmat(opt.local_weights, [1 mesh.K])*spdiag(mesh.cell_scale);
mesh.local_nodes = opt.local_nodes;
mesh.local_weights = opt.local_weights;

% Global indices of nodes that lie on faces
mesh.face_indices = [1 + (0:(mesh.K-1))*mesh.N; (1:mesh.K)*mesh.N];
mesh.face_indices = mesh.face_indices(:);

% Number each face, and assign an outward-pointing normal vector to it
mesh.face_normals = [-sign(diff(cells,1,2))].';
mesh.face_normals = [mesh.face_normals; -mesh.face_normals];
mesh.face_normals = mesh.face_normals(:);

vertices = mesh.cell_boundaries;

tol = 1e-12;
f2v = spalloc(mesh.K+1, 2*mesh.K, 2*mesh.K);
for q = 1:mesh.K;
    row = find(abs(mesh.cell_boundaries - mesh.cells(q,2))<tol);
    f2v(row, 2*(q-1)+2) = 1;
    row = find(abs(mesh.cell_boundaries - mesh.cells(q,1))<tol);
    f2v(row , 2*(q-1)+1) = 1;
end

f2f = (f2v.')*f2v - speye(2*mesh.K);
[a,b] = find(f2f);

% Each facial location has an "internal" and an "external" value. These indices
% mark the global node indices of the "internal" and "external" nodes. For the
% cells on the boundaries, the "external" node just points to the "internal"
% one; you can set boundary conditions manually later.

mesh.face_to_face = zeros([2*mesh.K 1]);
mesh.face_to_face(a) = b;
temp = find(mesh.face_to_face==0);
mesh.face_to_face(temp) = temp;

% Smallest dx:
mesh.dx = min(diff(mesh.local_nodes))*min(mesh.cell_scale);
