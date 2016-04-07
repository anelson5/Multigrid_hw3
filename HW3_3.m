%%%%%
%Homework 3 problem 3
%
%
%
%
close all
clear
clc
%Problem type
plots = 'on'; 

n = 2^7;
h = 1/n; %how many vcycles?
tol = 1e-10;
errorvect = zeros(1,1); 
%Pre and post smoothing
nu1 = 2;
nu2 = 1;

%functions
w_true = @(x) sin(pi*x).*(- x.^3 + x.^2); 
rhs = @(x,l) sin(pi*x).*(6*x - 2) - pi^2.*x.^2.*sin(pi*x).*(x - 1) +...
2*pi.*x.*cos(pi*x).*(3*x - 2) -l.*x.^2.*exp(-x.^2.*sin(pi*x).*(x - 1)).*sin(pi*x).*(x - 1);
%initial guess
w_0 = @(x) sin(pi*x).*(- x.^3 + x.^2); 

%Gamma
gam = 0; 

%Number of Newton Method steps
L = 3; 
%set up grid
x = 0:h:1; 
true= w_true(x); 
f = rhs(x,gam); 
w0 = w_0(x); 

%initial condition
w0 = 0.*x;

%set up interior
f_interior = f(2:end-1); 
w0_interior = w0(2:end-1); 

v = w0_interior'; 
%set up boundary vector
boundary = zeros(size(x)); 


error1 = w0 - true ;
err = max(abs(error1));
counter = 1; 
  errorvect(counter) = max(abs(w0-true)); 
err = 1; 
while err > tol
    [v,residual] = vcycle_nl(h,f_interior,v, nu1, nu2,gam,L,boundary);
    residual;
    err = abs(max(residual));
    %err = norm(residual)*h^(0.5);
    counter = counter+1; 
      errorvect(counter) = err; 
   
end
 plot(x,true); 
hold on; 
plot(x, [0 v 0],'rx-');





