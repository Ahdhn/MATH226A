function [normF,X] = Newton(Func, Jac, X0, MaxIt, LS)   
    %apply Newton method and return the function norm for each iteration 
    %along with the step 
    %@Func is the input function(s) for which we are seeking the roots    
    %@Jac is a function to evaluate the Jacobian inverse for Func    
    %@X0 is the initial guess     
    %@MaxIt is the max number of iterations     
    %@LS is a flag to use line search     
    
    x = X0;
    X = X0';
    normF=norm(double(Func(x)));
    n = 1;
    tol = 1e-6;%stopping tolerance 
    while n < MaxIt %don't go beyound iteration limit 
        xn = x- Jac(x)*Func(x);
        X = [X;xn'];
        while LS && norm(Func(xn)) > norm(Func(x))
             %if we are overshooting and Line Search 
             x = (x+xn)/2.0;
             xn = x- double(Jac(x))*double(Func(x));
        end
        normF(end+1) = norm(double(Func(xn)));               
        if abs(normF(end)) < tol || abs(norm(X(end)) - norm(X(end-1))) <tol
            %stopping criteria on the size of the norm and step size 
            break;
        end
        x = xn; %swap for next iteration  
        n=n+1;
    end
end