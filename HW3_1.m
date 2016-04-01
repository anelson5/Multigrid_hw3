%%Homework 3: Problem 1, solving -w`_xx - w_yy = rhs with two sets of boundary
%%conditions

n = 2^7; 
h = 1/n; 
v_cycle = 15; %how many vcycles?

%Pre and post smoothing
nu1 = 1; 
nu2 = 1; 

%Intial error vectors
error_max = zeros(v_cycle,1); 
error_2 = zeros(v_cycle,1); 
resid_max = zeros(v_cycle,1); 
resid_2 = zeros(v_cycle,1);

%Nodes, including boundary
[x,y] = ndgrid(0:h:1, 0:h:1); 
x
y
%Rhs based on exact solution, with negative right hand side

rhs= @(x,y) (25*pi^2.*exp((y - 1/2).^2 - (x - 1/2).^2).*cos(4*pi*y).*sin(3*pi*x)...
    - exp((y - 1/2).^2 - (x - 1/2).^2).*cos(4*pi*y).*sin(3*pi*x).*(2*y - 1).^2 - ...
    exp((y - 1/2).^2 - (x - 1/2).^2).*cos(4*pi*y).*sin(3*pi*x).*(2*x - 1).^2 +...
    6*pi*exp((y - 1/2).^2 - (x - 1/2).^2).*cos(3*pi*x).*cos(4*pi*y).*(2*x - 1)...
    + 8*pi*exp((y - 1/2).^2 - (x - 1/2).^2).*sin(3*pi*x).*sin(4*pi*y).*(2*y - 1)); 

w_true = @(x,y) sin(3*pi*x).*cos(4*pi*y).*exp(-(x-1/2).^2 + (y-1/2).^2);
true = w_true(x,y);
righthand = rhs(x,y);

% Boundary data
boundary = w_true(x,y); 
inside = zeros(size(boundary(2:end-1,2:end-1))); 
boundary(2:end-1,2:end-1) = inside;  
figure(99)
mesh(x,y,boundary)
figure(100)
mesh(x,y,true)
title('true solution')

w_0 = @(x,y)sin(3*pi*x).*cos(4*pi*y).*exp(-(x-1/2).^2 + (y-1/2).^2); 
figure(2)
mesh(x,y,w_0(x,y));
title('initial guess')


%Set up initial error vectors
diff = w_0(x,y)-w_true(x,y); 
error_max(1)=max(sum(abs(diff')));
norm2 = norm(diff,2).^2;
error_2(1) = (h^2*norm2)^0.5; 
w_interior = w_true(x,y); 
w_interior = w_interior(2:end-1,2:end-1); 

rhs = righthand(2:end-1, 2:end-1);
v = w_interior; 
true = true(2:end-1,2:end-1); 
x_plot = x(2:end-1,2:end-1);
y_plot = y(2:end-1,2:end-1);
true
boundary(2:end-1,2:end-1) = true
%u = RBGS(v, rhs, nu1, h,boundary); 
%igure(99)
%boundary(2:end-1,2:end-1) = u; 
%boundary
mesh(x,y,boundary)
for counter = 1:v_cycle
    diff = v-true; 
     error_max(counter)=max(sum(abs(diff')));
    norm2 = norm(diff,2).^2;
    error_2(counter) = (h^2*norm2)^0.5; 
    [v,residual] = vcycle2d(h,rhs,v, nu1, nu2,boundary);
    %figure(counter+50)
    % mesh(x_plot,y_plot, v)
end
figure(3)
mesh(x_plot,y_plot,residual)
title('residual'); 
figure(4)
mesh(x_plot,y_plot,true-v)
title('error'); 
figure(6)
boundary(2:end-1,2:end-1) = v; 
mesh(x,y,boundary)
title('output from multigrid'); 
xlabel('x'); 
ylabel('y'); 