function [output] = nonlinear_op(v,h,boundary,gam)
%Nonlinear operatotr
%v - input into operator, interior
%boundary - boundary vector
%gam- gamma assosiated with operator w'' + gamma*w*exp(x) = f
n = length(v);
boundary(2:end-1) = v;
v_wb = boundary;
Lv = zeros(1,n+2); 
for j = 2:n+1
    Lv(j) = 1/h^2*(-v_wb(j-1) + 2*v_wb(j) - v_wb(j+1)) + gam*v_wb(j)*exp(v_wb(j)); 
end

output = Lv(2:end-1); 


end

