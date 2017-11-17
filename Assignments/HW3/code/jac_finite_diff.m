%https://www.mathworks.com/matlabcentral/answers/28066-numerical-jacobian-in-matlab
clc;
clear;
close all;

% Jacobian functor
J = @(x,h,F)(F(repmat(x,size(x'))+diag(h))-F(repmat(x,size(x'))))./h';
% Your function
f = @(x)[x(1)^2 + x(2)^2; x(1)^3];
% Point at which to estimate it
%x = [50;13];
x = [5; 2];

% Step to take on each dimension (has to be small enough for precision)
h = 1e-5*ones(size(x));

% Compute the jacobian
%J(x,h,f)

global x11s x12s f11s f12s; 
syms x11s x12s f11s f12s;
f11s = (5-2*x11s)/(2*x12s-3);%function in symbolic 
f12s = (5-2*x12s)/(2*x11s-3);%to get the jacobian 

estJac(x,h, f)



function J = estJac(x,h,F)

    F1 = F(repmat(x,size(x'))+diag(h));
    
    F2 = F(repmat(x,size(x')));    
    
    J = F1 - F2;
    
    J = J./h';
    
end

function val = Func1(x)
    global x11s x12s f11s f12s
    %evluate the function on x
    %x is a vector
    x11s = x(1);
    x12s = x(2);
    val(1) = subs(f11s);
    val(2) = subs(f12s);
    val = double(val)';
end