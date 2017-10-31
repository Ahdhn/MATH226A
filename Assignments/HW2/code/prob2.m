clc;
clear;
close all;
disp('Math 226A - HW2 - Problem #2:');
start = 1;
max_size = 250;

i=1;
for m_size=start:max_size
    [error(i), cond_num(i)] = trial(m_size);
    i=i+1;
end

scatter([start:max_size],error);
set(gca,'yscale','log');
ylabel('Relative Error','FontSize',15);
xlabel('Matrix Size','FontSize',15);



hold on
yyaxis right
ylabel('Condition Number','FontSize',15);

scatter([start:max_size],cond_num,'*');
set(gca,'yscale','log');
set(gca,'xscale','log');
function [relative_error, condition_num] = trial(m)
    %generate an m x m random square matrix A ,
    %generate a m-length random vector Z,
    %compute b = Az
    %solve Ax = b for x
    %return the relative error |x-z|/|z| and condition number of A
    A = randn(m);
    z = randn(m,1);
    b = A*z;
    x = A\b;
   condition_num = cond(A);
   relative_error = norm(x-z)/norm(z);
end