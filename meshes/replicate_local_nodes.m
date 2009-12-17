function[global_nodes] = replicate_local_nodes(local_nodes, cells)
% replicate_local_nodes -- creates global set of nodes from local node template
%
% global_nodes = replicate_local_nodes(local_nodes, cells)
%
%     Given a size K x 2 vector cells defining K intervals, and a length-N
%     vector local_nodes, this function returns an N x K matrix where column k
%     corresponds to the vector of nodes local_nodes mapped from [-1,1] ----->
%     [cells(k,1) cells(k,2)].

persistent spdiag
if isempty(spdiag)
  from labtools import spdiag
end

local_nodes = local_nodes(:);
N = length(local_nodes);
K = size(cells,1);

scales = diff(cells,1,2)/2;
shifts = mean(cells,2);

global_nodes = local_nodes*ones([1 K]);
global_nodes = global_nodes*spdiag(scales) + ones([N 1])*shifts.';
