clc;
clear;
close all;
disp('Math 226A - HW1 - Problem #4:');

len = 20;
my_coef = poly(1:len);
myfunc = @(x)(dot(my_coef,fliplr(x.^(0:len))));
for n=1:20
    myroots(n)= fzero(myfunc,n);
end


