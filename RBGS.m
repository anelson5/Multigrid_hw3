function [v_new,residual] = RBGS(v, rhs, nu, h,boundary)
%Relax - Performs weighted red black Gauss Seidel relaxation for the 2D
%problem with -w_xx - w_yy = f(x,y) on the unit square
% v -- V cycle output (in matrix for 2D problem), 
% rhs -- right hand side, as a vector f_{j,l} on the 2D grid interior nodes
% h -- grid stepsize
% nu -- number of relaxations 

%change v into a matrix, v is (n-1)^2 we want v to be (n-1)x(n-1)
n = size(v);
n = n(1); 
[X,Y] = ndgrid(0:h:1,0:h:1);

%Boundary data; 
boundary_matrix = boundary; 
boundary_matrix(2:n+1,2:n+1) = v; 
v_wb = boundary_matrix;
f = rhs; 
for k = 1:nu

 
     %Odd nodes
    for l = 2:n+1
        for j = 2:n+1                    
          if mod(j+l,2) == 1
          v_wb(j,l) = 1/4*(v_wb(j,l-1) + v_wb(j-1,l) + v_wb(j,l+1) + v_wb(j+1,l) + h^2*f(j-1,l-1));
          end
        end
    end
      %Even nodes
    for l = 2:n+1
        for j = 2:n+1
            if mod(j+l,2) == 0
           v_wb(j,l) = 1/4*(v_wb(j,l-1) + v_wb(j-1,l) + v_wb(j,l+1) + v_wb(j+1,l) + h^2*f(j-1,l-1));
            end
        end
    end
end
v_new = v_wb(2:n+1, 2:n+1);
 
Lv_new = zeros(n+2,n+2); 
for l = 2:n+1
    for j = 2:n+1
        Lv_new(j,l) = (v_wb(j-1,l) - 4*v_wb(j,l) + v_wb(j+1,l) + ... 
            + v_wb(j,l-1) + v_wb(j,l+1))/h^2; 
    end
end
Lv = Lv_new(2:n+1, 2:n+1);
residual_matrix = f+Lv;

residual = residual_matrix;
%mesh(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),residual)
