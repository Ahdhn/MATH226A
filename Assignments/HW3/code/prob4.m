%HW3 - problem 3
clc;
clear;
close all;
disp('Math 226A - HW3 - Problem #4:');
%A = [1 1 4; 1 4 2; 1 4 2; 1 1 0];
%[Q1,R1] = GS(A);
%[Q2,R2] = MGS(A);
num_points = 50;
poly_deg = 11;

%sample the function 
sams = linspace(0,1,num_points);
func = cos(4.*sams);


%computes all the powers once
p = zeros(2*poly_deg +1,1);
for k=1:num_points
    for n=1:2*poly_deg +1
        p(n) = p(n) + sams(k)^(n-1);
    end
end
%construct the normal matrix
A = zeros(poly_deg +1,poly_deg+1);
for k = 1:poly_deg +1
    for n = 1:poly_deg +1
        A(k,n) = p(k+n-1);
    end
end
%construct the right head side b
b = zeros(poly_deg+1,1);
for k=1:num_points
    for n = 1:poly_deg +1
        b(n) = b(n) + func(k)*sams(k)^(n-1);
    end
end

%solve using classic GS
[Q_GS, R_GS]=GS(A);
c_GS = inv(R_GS)*transpose(Q_GS)*b;
norm_residual_GS = norm(func - polyval(c_GS,sams)); 
error_GS = GetError(c_GS);


%solve using modified GS
[Q_MGS, R_MGS] = MGS(A);
c_MGS = inv(R_MGS)*transpose(Q_MGS)*b;
norm_residual_MGS = norm(func - polyval(c_MGS,sams)); 
error_MGS = GetError(c_MGS);

%solve using Housholder 
[Q_H, R_H]=qr(A);
c_H = inv(R_H)*transpose(Q_H)*b;
norm_residual_H = norm(func - polyval(c_H,sams)); 
error_H = GetError(c_H);

function err = GetError(my_c)
    err(1) = abs(my_c(1) - 1.000000000996606);
    err(2) = abs(my_c(2) - -4.227430949815150e-07);
    err(3) = abs(my_c(3) - -7.999981235683346e+00);
    err(4) = abs(my_c(4) - -3.187632625738558e-04);
    err(5) = abs(my_c(5) - 1.066943079610163e+01);
    err(6) = abs(my_c(6) - -1.382028878048870e-02);
    err(7) = abs(my_c(7) - -5.647075625417684e+00);
    err(8) = abs(my_c(8) - -7.531602738192263e-02);
    err(9) = abs(my_c(9) - 1.693606966623459e+00);
    err(10) = abs(my_c(10) - 6.032106743884792e-03);
    err(11) = abs(my_c(11) - -3.742417027133638e-01);
    err(12) = abs(my_c(12) - 8.804057595513443e-02);
end
%solve using matlab's least square solver
%c_matlab = A\b;
%norm_residual_matlab = norm(func - polyval(c_matlab,sams)); 


