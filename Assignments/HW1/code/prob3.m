clc;
clear;
close all;
disp('Math 226A - HW1 - Problem #3:');



dt1 = 0.1;
t1 = 0;
Nsteps1 = 864000;
for j=1:Nsteps1
    t1=t1+dt1;
    t1_exact = j*dt1;
end


dt2 = 0.125;
t2 = 0;
Nsteps2 = 691200;
for i=1:Nsteps2
    t2=t2+dt2;
    t2_exact = i*dt2;
end