clc;
clear;
close all;
disp('Math 226A - HW1 - Problem #1:');

disp('a) For n=10, relative error are:');
n = 10;  
disp_error(0.5, ApproxExpFunc(0.5,n));
disp_error(-0.5, ApproxExpFunc(-0.5,n));
disp_error(30*pi, ApproxExpFunc(30*pi,n));
disp_error(-30*pi, ApproxExpFunc(-30*pi,n));
disp('For n=40, relative error are:');
n = 40;     
disp_error(0.5, ApproxExpFunc(0.5,n));
disp_error(-0.5, ApproxExpFunc(-0.5,n));
disp_error(30*pi, ApproxExpFunc(30*pi,n));
disp_error(-30*pi, ApproxExpFunc(-30*pi,n));


disp('b) For n=10, relative error are:');
n = 10;  
disp_error(0.5, AccurateExpFunc(0.5,n));
disp_error(-0.5, AccurateExpFunc(-0.5,n));
disp_error(30*pi, AccurateExpFunc(30*pi,n));
disp_error(-30*pi, AccurateExpFunc(-30*pi,n));

function accurate_exp = AccurateExpFunc(x,n)
    %More accurate approximation for the exp function
    %it exploites the fact that exp(x) = exp(x/2+x/2) = exp(x/2)*exp(x/2)
    %with in mind that ApproxExpFunc() give good approximation 
    %(relative err < 10^-13) for x<0.2
    
    if abs(x) < 0.2
       accurate_exp = ApproxExpFunc(x, n);
    else 
        accurate_exp = AccurateExpFunc(x/2, n);
        accurate_exp = accurate_exp*accurate_exp;
    end
end
function approx_exp = ApproxExpFunc(x, n)
    %Approximate the exponential using truncated series    
    approx_exp = 1;  %first term is 1 (x^(0)/0!)
    j_factorial = 1; %start with factorial of 0
    for j=1:n
        j_factorial = j_factorial*(j);
        approx_exp = approx_exp + (x^(j))/j_factorial;    
    end
end
function disp_error(x,approx_exp)
    %display the relative error between exp(x) and approx_exp
    fprintf('x=%f   --> %e\n', x ,abs(exp(x) - approx_exp)/exp(x));
end

