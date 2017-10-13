clc;
clear;
close all;
disp('Math 226A - HW1 - Problem #4:');

len = 20;
my_coef = poly(0:len);
myfunc = @(x)(dot(my_coef,fliplr(x.^(0:len))));

%i=1;
%for n=1:.01:50
 %   y(i) = myfunc(n);
 %   i=i+1;
%end
%plot(1:i-1, y);

for n=1:20
    myroots(n)= fzero(myfunc,n);
end
