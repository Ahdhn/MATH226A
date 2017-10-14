clc;
clear;
close all;
disp('Math 226A - HW1 - Problem #2:');


values = zeros(2,10);
for i=1:12
    x = power(10,i-1);
    values(1,i) = atan(x+1) - atan(x);
   % my_t1= log((x+1-1i)/(x+1+1i));
   % my_t2 = log((x-1i)/(x+1i));
   % my_tt = my_t1 - my_t2;
   % values(2,i) = (1/2i)*(my_tt);
   values(2,i) = (1/2i)*(log((x+1-1i)/(x+1+1i)) /log((x-1i)/(x+1i)));
end


syms s;
f = 5*s/s;