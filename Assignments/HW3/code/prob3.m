%HW3 - problem 3
clc;
clear;
close all;
disp('Math 226A - HW3 - Problem #3:');

%%%%%%%%%%%%%%%%%%% Part a)
global x11s x12s j1s f11s f12s; 
syms x11s x12s j1s f11s f12s;
f11s = (5-2*x11s)/(2*x12s-3);%function in symbolic 
f12s = (5-2*x12s)/(2*x11s-3);%to get the jacobian 
j1s = [jacobian(f11s,[x11s, x12s]);jacobian(f12s,[x11s, x12s])];
j1s = inv(j1s);%the inverse
disp('Results for part a) x0=(1 1), no Line Search');
[normF, X] = Newton(@Func1, @myJac1, [1 1]', 3 ,true);  %initial guess (1,1), no LS


v=[.25,.5:2:40];
xr=-5:0.2:5;
n=length(xr);
z=zeros(n,n);
for i=1:n
    for j=1:n        
        z(i,j)=norm(Func1([xr(i),xr(j)]));
    end
end
figure(1);
contourf(xr,xr,z,v);
hold 
plot(X(:,1),X(:,2),'-*');

function j = myJac1(x)
    global x11s x12s j1s
    x11s = x(1);
    x12s = x(2);   
    norm(double(subs(j1s)));%get an idea of the matrix singulairty
    j = double(subs(j1s));
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