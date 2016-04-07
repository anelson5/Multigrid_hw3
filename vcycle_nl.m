function [v,error,residual] = vcycle_nl(h,f,v_0, nu1, nu2,gam,L,boundary)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
% h - step size
% f_array - cell array of rhs
% v_arry - cell array of v
% nui - relaxation parameter. 

if length(v_0) ==1
%if levels == 1 
    %direct solve with nonlinear Newton's method
    v = nm(v_0,f,gam,boundary,L,h);
    error = v-v_0; 
    residual = 0; 
else
    [v_nu1, residual] = nonlinear_relax(v_0, f, nu1, h,boundary,gam);
    v= v_nu1; 
    v_coarse = restrict(v);
    r_coarse = restrict(residual); 
    boundary_coarse = zeros(length(v_coarse)+2,1);
    f_coarse = r_coarse' + nonlinear_op(v_coarse,2*h,boundary_coarse,gam);
    [v_cycle,error_coarse,~] = vcycle_nl(2*h,f_coarse, v_coarse, nu1, nu2,gam,L,boundary_coarse);
    if iscolumn(v_cycle) == 0
        v_cycle = v_cycle';
    end
     
    if iscolumn(v_0) == 0
        v_0 = v_0';
    end
     if iscolumn(v_nu1) == 0
        v_nu1 = v_nu1';
     end
    
     
    
    %v_0 or v_nu1 + interpolate(error)?
    v = v_nu1 + interpolate(error_coarse); 
    
    [v_nu2, residual] = nonlinear_relax(v, f, nu2, h,boundary,gam);
    v= v_nu2; 
    if iscolumn(v_nu2) == 0
        v_nu2 = v_nu2'; 
    end
    
    if iscolumn(v_nu1) == 0
        v_nu1 = v_nu1'; 
    end
    
    %is this v_nu2-v_nu1? or v_nu2 - v_0?
 
    error = v_nu2 - v_nu1;
    
 
end

