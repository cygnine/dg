function[um] = u_mean(u_minus, u_plus)
% u_mean -- computes the mean of two inputs
%
% um = u_mean(u_minus, u_plus)
%
%     Just a fancy wrapper for an average of two like-sized arrays. 

um = 1/2*(u_minus + u_plus);
