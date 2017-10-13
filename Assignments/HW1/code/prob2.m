clc;
clear;
close all;
disp('Math 226A - HW1 - Problem #2:');

x = 10;
values = zeros(2,10);
i=1;
while x<10^11
    values(1,i) = atan(x+1) - atan(x);
    values(2,i) = (1/2i)*(log((x+1-1i)/(x+1+1i)) - log((x-1i)/(x+1i)));
    x = x*10;
    i = i+1;
end

