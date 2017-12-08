%HW4 - Problem 2 

clc;
clear;
close all;

%p = [21 41 81 161 321 641 1281];
%err_nat = [0.0031827 2.7797e-04 1.6107e-05 1.6142e-06 4.0366e-07 1.009220670436517e-07 2.523093434875223e-08];
%err_knot = [0.0031827 2.7797e-04 1.6107e-05  9.6744e-07 5.9822e-08 3.728698017013699e-09 2.328839343590516e-10];
%plot (p,err_nat, p, err_knot,'LineWidth',2);
%legend('Natural Spline', 'Not-a-knot Spline');
%set(gca,'yscale','log');
%set(gca,'xscale','log');
%xlabel('Number of points');
%ylabel('Error');


x0 = -1;
x1 = 1;
spacing = 0.0015625;

%myFun = @(x)(cos(2*pi*x));
myFun = @(x)( 1./(1 + 25.*x.^2));

x = (x0:spacing:x1)';
y = myFun(x);

P_not_a_knot = naturalspline(x,y,true);
P_natural = naturalspline(x,y,false);

z = (x0 : spacing/100 : x1)';
z_natural  = zeros(length(z),1);
z_not_a_knot = zeros(length(z),1);

for l=1:1:length(z)
    z_natural(l) = Eval(P_natural,x,z(l));
    z_not_a_knot(l) = Eval(P_not_a_knot,x,z(l));    
end

x_true = (x0:spacing/100:x1)';
y_true = myFun(x_true);

plot(x_true,y_true,'LineWidth',2);
hold on
plot(z,z_natural,'o');
plot(z,z_not_a_knot,'*');
legend('Ground Truth','Natural','Not-a-knot');
xlabel('X');
ylabel('F(x)');

if  length(z_natural) == length(y_true)
    error_natural = max(abs(z_natural-y_true))
end
if  length(z_not_a_knot) == length(y_true)
    error_not_a_knot = max(abs(z_not_a_knot-y_true))
end

function val = Eval(P, x, z)
    %P is the coefficients of the natural cubic spline 
    %x defines the intervals through which P was constructed 
    %z is a scalar within x to evelate S(z)
    
    if z>max(x) || z<min(x)
        print('z is outside x');
    end
    
    %find in which interval z is located 
    for l= 1:1:length(x)-1 
        if x(l)<=z && x(l+1)>=z
            break;
        end
    end
    val = P(l,1) + P(l,2)*(z-x(l)) + P(l,3)*(z-x(l))^2 + P(l,4)*(z-x(l))^3;
    
end