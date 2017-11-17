%HW3 - problem 3
clc;
clear;
close all;
disp('Math 226A - HW3 - Problem #3:');

global my_factor my_g

N = 100;%number of mesh points 
delta_x = 1/(N+1); %spacing 
my_g = 0.25; 
my_eps = 1e-3;
my_factor = my_eps / (delta_x*delta_x); %matrix multiplier 
x = 1:N; %mesh points 

%first initial guess is 0.25 for all points 
init1 = my_g*ones(length(x)-2,1);%first initial conditions 

%second initial guess is 0.5*x^2*(my_g/my_eps)
init2=zeros(length(x)-2,1);%second initial conditions 
for n = 2:N-1
    init2(n-1) = -0.5*(my_g/my_eps)*(n)^2;
end

[normF1, myU1] = Newton(@Func, @myJac, init1, 1000 ,true); 
figure 
plot(myU1(end,:));
xlabel('x'); ylabel('u');
title('Solution with first initial guess');

[normF2, myU2] = Newton(@Func, @myJac, init2, 1000 ,true); 
figure 
plot(myU2(end,:));
xlabel('x'); ylabel('u');
title('Solution with second initial guess');


function val = Func(U)
    %evaluate the discretized problem at U
    global my_factor my_g    
    val = zeros(length(U) , 1); %the solution should be a column     
   for n=1:length(U)
       
       %u_{j-1}
       if n == 1
           u_mins = 0;
       else 
           u_mins = U(n-1);
       end
       
       %u_{j+1}
       if n == length(U)
           u_plus = 0;
       else 
           u_plus = U(n+1);
       end
       
       %u
       u = U(n);       
       val(n) = my_factor*(u_mins -2.0*u + u_plus)...
           - exp(-abs(u))*u + my_g;
   end
end
function jac = myJac(U)
    %compute the tridiagonal jacobian using finite difference
    %as mentioned in Kelley book 
    f0 = Func(U);
    n = length(U);
    jac = zeros(n,n);
    dv = zeros(n,1);
    epsnew = 1.e-7;
    for ip = 1:n
        delr(ip) = min([1+1+ip,n]);
        ih(ip)=min([ip+1,n]);
        il(ip)= max([ip-1,1]);
    end
    for is = 1:delr(1)
        ist = is;
        pt = zeros(n,1);
        while ist<= n
            pt(ist) =1;
            ist = delr(ist)+1;
        end
        U1 = U+epsnew*pt;
        f1 = Func(U1);
        dv = (f1-f0)/epsnew;
        ist = is;
        while ist <= n
            ilt = il(ist);
            iht = ih(ist);           
            jac(ilt:iht,ist) = dv(ilt:iht);
            ist = delr(ist)+1;
        end
    end
    jac = inv(jac);
end
