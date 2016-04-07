function [v_new, residual] = nonlinear_relax(v, rhs, nu, h,boundary, gam)
%NONLINEAR_RELAX: Performs nonlinear GS relaxation for the Full
%approximation scheme problem, number 3 Will use Newtons method
%
%
n = length(v);
N = length(boundary);
f = rhs;

%Number of Newton Method iterations
L = 3;
%

for k = 1:nu
    v_new = nm(v,f,gam,boundary,L,h);
    
end

applied = nonlinear_op(v,h,boundary,gam);

if iscolumn(applied) == 0
    applied = applied';
end
if iscolumn(f) == 0
    f = f';
end

residual = f - applied; 


end

