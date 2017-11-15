function normF = Newton(Func, Jac, X0, MaxIt, LS)   
    %apply Newton method and return the function norm for each iteration 
    %Func is the input function(s) for which we are seeking the roots    
    %Jac is a function to evaluate the Jacobian for Func    
    %X0 is the initial guess     
    %MaxIt is the max number of iterations     
    %LS is a flag to use line search     
    x = X0;
    normF=norm(double(Func(x)));
    n = 1;
    while n < MaxIt         
        xn = x- Jac(x)*Func(x)
        while LS && norm(Func(xn)) > norm(Func(x))
             %if we are overshooting and line search 
             x = (x+xn)/2.0;
             xn = x- double(Jac(x))*double(Func(x));
        end
        normF(end+1) = norm(double(Func(xn)))
        x = xn; %swap for next iteration  
        n=n+1;
    end
end