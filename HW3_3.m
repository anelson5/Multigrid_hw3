%%%%%
%Homework 3 problem 3
%
%
%
%

rhs = @(x,l) sin(pi*x).*(6*x - 2) - pi^2.*x^2.*sin(pi*x).*(x - 1) +...
2*pi*x.*cos(pi*x).*(3*x - 2) - l.*x.^2.*exp(-x^.2.*sin(pi*x).*(x - 1)).*sin(pi*x).*(x - 1)