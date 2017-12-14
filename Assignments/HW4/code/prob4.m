%HW4 - Problem 4
clc;
clear;
close all;

%myFun = @(x)(1/(x+0.01));

[q,new_count]=adapt(@myFunc,0,1, 1e-10);

[Q, Fcnt] = quad(@myFunc, 0 , 1 , 1e-10);

itr =1;
truth = log(1.01) - log(0.01);
for n=4:2:100000
    error(itr) = abs(Simpson(@myFunc, n, 0,1) -truth);
    if error(itr) < 10^-10
        disp('Min Error is=');
        error(itr)
        disp ('Number of points is');
        n
        break;
    end    
    itr = itr+1;  
end


function result = Simpson(Func, n, a,b)
    %apply composite Simpson's rule to approximate the intergation 
    %of Func between a and b using n-subintervals 
    result=Func(a);
    h = (b-a)/n;
    for int =1:1:n-1
        if mod(int,2) == 0
            factor = 2;%even
        else
            factor = 4; %odd
        end
        result = result + factor * Func(a + int*h);
    end
    result = result + Func(b);
    result = result * h/3;
end


function y = myFunc(x)
    y = 1./(x+0.01);
end