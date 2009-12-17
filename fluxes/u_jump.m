function[uj] = u_jump(u_minus, u_plus, normal_minus, normal_plus)
% u_jump -- Computes the jump [[u]] of a multivalued function
%
% uj = u_jump(u_minus, u_plus, normal_minus, normal_plus)
%
%     The inputs u_minus and u_plus are like-sized column vectors containing
%     function values on the interior ("minus") and exterior ("plus") of an
%     interface. The normal_minus and normal_plus arrays are N x D arrays, where
%     N is the length of u_minus, and D is the spatial dimension of the
%     approximation. Each row of normal_minus represents the D-dimensional
%     normal vector at the interface. Minus is inward-pointing, plus is
%     outward-pointing. 
% 
%     This implements the "scalar" version of the jump operator (u_minus is a
%     scalar), and thus returns a D-dimensional output, uj is N x D.

D = size(normal_minus,2);
N = length(u_minus);

uj = zeros([N D]);
uj = (u_minus*ones([1 D])).*normal_minus + (u_plus*ones([1 D])).*normal_plus;
