%HW4 - Problem 3
clc;
clear;
close all;

myFun = @(x)(exp(x));
truth = exp(1.0) - exp(0.0);
attempts = 5;
error = zeros(attempts,3);
n=2;
for itr=1:1:attempts 
    error(itr,1) = abs(Trapezoidal(myFun, n, 0,1) - truth);
    error(itr,2) = abs(Simpson(myFun, n, 0,1) - truth);
    error(itr,3) = abs(GaussianQuad(myFun, n, 0,1)-truth);
    n =n*2;    
end

plot(error(:,1),'LineWidth',2);
hold on
plot(error(:,2),'LineWidth',2);
plot(error(:,3),'LineWidth',2);

legend('Trapezoidal Rule', 'Simpson''s Rule', '3-pnt Gaussian quad');
set(gca,'yscale','log');
set(gca,'xscale','log');
xlabel('Number of points');
ylabel('Error');

function result = Trapezoidal(Func, n, a,b)
    %apply composite trapezoidal rule to approximate the intergation 
    %of Func between a and b using n-subintervals 
    result = Func(a)/2.0;
    h = (b-a)/n;
    for int = 1:1:n-1
        result = result + Func(a + int*h);
    end
    result = result + Func(b)/2.0;
    result =result * h;
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


function result = GaussianQuad(Func, n, a,b)
    %apply Gaussian Quadrature rule to approximate the intergation 
    %of Func between a and b using n-subintervals 
    result=0;    
    h = (b-a)/n;
    for int =0:1:n-1
       %for each interval apply the 3-point GQ
       myA = a + int*h;
       myB = myA + h;
       
       myH = (myB-myA)/18;
       
       p0 = ((myA + myB)/2.0) - sqrt(3.0/5.0)*((myB- myA)/2.0);
       p1 = ((myA + myB)/2.0);
       p2 = ((myA + myB)/2.0) + sqrt(3.0/5.0)*((myB- myA)/2.0);
       
       result = result + myH*(5.0*Func(p0) + 8.0*Func(p1) + 5.0*Func(p2));       
    end
end