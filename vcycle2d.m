function [v,residual] = vcycle2d(h,f,v, nu1, nu2, boundary)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
% l - level numbers
% h - step size
% f_array - cell array of rhs
% v_arry - cell array of v
% nui - relaxation parameter.

[v_new, residual] = RBGS(v, f, nu1, h,boundary);
v = v_new;
[~,n] = size(v);
if n == 1
    %if levels == 1
    %direct solve
    
    v = directsolve2d(v,f,h,boundary);
    %go to next relax
    
else
    %figure(100)
   % mesh(x1(2:end-1,2:end-1),y1(2:end-1,2:end-1),residual); 
    f_coarse = restrict2d(residual);
    %figure(101)
%    mesh(x(2:end-1,2:end-1),y(2:end-1,2:end-1),f_coarse); 
    vcoarse = zeros(size(f_coarse));
     boundary_coarse = zeros(size(boundary(1:2:end,1:2:end)));
    [v_cycle] = vcycle2d(2*h,f_coarse,vcoarse, nu1, nu2,boundary_coarse);
end

if n ~= 1
    v = v + interpolate2d(v_cycle,boundary_coarse);
end
[v_new, residual] = RBGS(v, f, nu2, h, boundary);
v= v_new;

end


