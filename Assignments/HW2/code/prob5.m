clc;
clear;
close all;
disp('Math 226A - HW2 - Problem #5 (a):');
%Solving the cyclic tri-diagonal system to equations 

%Some test values 
n = 5; %Num row/colums
b = [-2 -2 -2 -2 -2]; %Diagonal
a = [2 1 1 1 1]; %Sub-diagonal (added upper right corner as first entry)
c = [1 1 1 1 2]; %Upper-diagonal (added lower left corner as last entry)
f = [1 1 1 1 1];%RHS

alpha = c(n); %Lower left corner
beta = a(1); %Right upper corner
bb = b; %Modified diagonal vector
gamma = -b(1);
bb(1) = b(1) - gamma;
bb(n) = b(n) - alpha*beta/gamma;
x = myTridiagonalSolver(n,bb,a,c,f); %solve  A.x = f for u 
u(1) = gamma;
u(n) = alpha;
z = myTridiagonalSolver(n,bb,a,c,u); %solve A.z= x for z
fact = (x(1) + beta *x(n)/gamma)/(1.0 +z(1) + beta*z(n)/gamma);
for i=1:n
    x(i) = x(i) - fact*z(i);
end


disp('The solution of X =');
disp(x);

function x = myTridiagonalSolver(n,b,a,c,f)
    x = zeros(n,1);
    w = zeros(n,1);
    d =b(1);
    x(1) = f(1)/d;
    for i = 2:n
         w(i-1) = c(i-1)/d;
         d = b(i) - a(i)*w(i-1);
         x(i) = (f(i)-a(i)*x(i-1))/d;
    end
    for i = n-1:-1:1
        x(i) = x(i) - w(i)*x(i+1);    
    end
end
