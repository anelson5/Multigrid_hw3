function [v_new] = interpolate2d(v_old,boundary)
%INTERPOLATE - interpolates fine grid nodes from the coarse grid
%   Set up interpolation matrix 
%   Set up matrix of (n-1)x(n/2-1)
%   Form (n/2-1)x(n/2-1) and n/2x(n/2-1) 

%Boundary data; 
n = size(v_old);
n = n(1);
boundary_matrix = boundary; 

boundary_matrix(2:n+1,2:n+1) = v_old;

v_wb = boundary_matrix;
[~,n1] = size(boundary_matrix); 
N = 2*n+2;
%Make the new finer grid 
v_new(2:2:N-2,2:2:N-2) = v_wb(2:N/2,2:N/2); %even j, even l
v_new(1:2:N-1,2:2:N-2) = 1/2*(v_wb(2:N/2+1,2:N/2) + v_wb(1:N/2,2:N/2)) ;
v_new(2:2:N-2,1:2:N-1) = 1/2*(v_wb(2:N/2,2:N/2+1) + v_wb(2:N/2,1:N/2)) ;
v_new(1:2:N-1,1:2:N-1) = 1/4*( v_wb(2:N/2+1,2:N/2+1) + v_wb(2:N/2+1,1:N/2)+...
    + v_wb(1:N/2,2:N/2+1) + v_wb(1:N/2,1:N/2));

end

