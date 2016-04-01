function [v_new] = directsolve2d(v,rhs,h,boundary)
%Direct solve Summary of this function goes here
%   Detailed explanation goes here
[~,n] = size(v); 
L = boundary;

v_new = 1/4*(L(1,2) + L(2,1) + L(3,2) + L(2,3) + h^2*rhs); 

end

