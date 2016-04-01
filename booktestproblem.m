%%Book test problem solving -w_xx - w_yy = rhs with two sets of boundary
%%conditions

n = 2^6; 
h = 1/n; 
v_cycle = 1; %how many vcycles?

%Pre and post smoothing
nu1 = 2; 
nu2 = 1; 

%Initial error vectors
error_max = zeros(v_cycle,1); 
error_2 = zeros(v_cycle,1); 
resid_max = zeros(v_cycle,1); 
resid_2 = zeros(v_cycle,1);

%Nodes, including boundary
[x,y] = ndgrid(0:h:1, 0:h:1); 

%Rhs based on exact solution, witH POS right hand side
rhs = @(x,y) -2*((1-6*x.^2).*y.^2.*(1-y.^2) + (1-6*y.^2).*x.^2.*(1-x.^2)); 
%rhs = @(x,y) 0*x; 
w_true = @(x,y) (x.^2-x.^4).*(y.^4-y.^2);
%w_true = @(x,y) 0.*x + 1; 

true = w_true(y,x);
righthand = rhs(y,x);


% Boundary data
x_ = linspace(0,1,n+1);
top_b = zeros(n+1,1); 
bottom_b = top_b; 
left_b = zeros(1,n+1);
right_b = left_b; 

boundary = zeros(n+1,n+1); 
boundary(1,:) = top_b; 
boundary(end,:)=bottom_b; 
boundary(:,1) = left_b; 
boundary(:,end) = right_b; 

%Look at boundary matrix, does it make sense
figure(1)
mesh(y,x,boundary)
xlabel('x')
ylabel('y')

%Initial guess that satisfies the boundary conditions
w_0 = @(x,y) 0.*x + 0;
figure(2)
mesh(x,y,w_0(x,y));
xlabel('x')
ylabel('y')

%Set up initial error vectors
diff = w_0(x,y)-w_true(x,y); 
error_max(1)=max(sum(abs(diff')));
norm2 = norm(diff,2).^2;
error_2(1) = (h^2*norm2)^0.5; 
w_interior = w_0(x,y); 
w_interior = w_interior(2:end-1,2:end-1); 


rhs = rhs(y,x); 
rhs = rhs(2:end-1, 2:end-1);
v = w_interior; 
true = true(2:end-1,2:end-1); 
for counter = 1:v_cycle
    diff = v-true; 
     error_max(counter)=max(sum(abs(diff')));
    norm2 = norm(diff,2).^2;
    error_2(counter) = (h^2*norm2)^0.5; 
    [v,residual] = vcycle2d(h,rhs,v, nu1, nu2,boundary);
    
    resid_max(counter)=max(sum(abs(residual')));
    norm2 = norm(residual,2).^2;
    resid_2(counter) = (h^2*norm2)^0.5; 
   
end
figure(3)
x_plot = x(2:end-1,2:end-1);
y_plot = y(2:end-1,2:end-1);
mesh(x_plot,y_plot,residual)
title('residual'); 
figure(4)
mesh(x_plot,y_plot,true+v)
title('error'); 
figure(5)
mesh(x_plot,y_plot, -v)

figure(6)
mesh(x_plot,y_plot,true); 