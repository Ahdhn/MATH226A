%HW4 - Problem 4
clc;
clear;
close all;

myFun = @(x)(1/(x+0.01));

[q,new_count]=adapt(myFun,0,1, 1e-10);

[Q, Fcnt] = quad(myFun, 0 , 1 , 1e-10);