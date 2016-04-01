function [v_new] = restrict2d(v_old)
%RESTRICT - performs a restriction operation from the fine grid to the
%coarse grid. The type is full weighting. 
% 
% v_old is row by row: v_old = v11 v12 v13, ..., v21, v22, v23, ... 
% v_old is (n-1)^2 x 1, v_new is (n/2-1)^2 x 1

[m,n] = size(v_old);
N_coarse = (n+1)/2-1;
v = v_old; 
v_new = zeros(N_coarse, N_coarse); 
for l = 1:N_coarse
    for j = 1:N_coarse
        v_new(j,l) = 1/4*v(2*j,2*l) + 1/8*(v(2*j,2*l+1) + v(2*j,2*l-1) + v(2*j-1,2*l) + ... 
            + v(2*j+1,2*l)) + 1/16*(v(2*j-1,2*l-1) + v(2*j-1,2*l+1) + v(2*j+1,2*l-1) ...
        + v(2*j+1, 2*l+1));
    end
end


end