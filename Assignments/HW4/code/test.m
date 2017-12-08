clc;
clear;

n= 10;
X = 0;
for k=1:1:n-1
   X(k) = (2*k-1) /(2*(n-1));
end
Y = 0;
for k=1:1:n+1
   Y(k) = (2*k-1) /(2*(n+1));
end

xcos = cos(pi.*X);
ycos = cos(pi.*Y);

plot(xcos)

res =1.0;
n= 50;
K=0;
for J=1:n
    res = res*(sin((pi*J)/(2*n)))^2;
end 

res = 1/res;

for J=0:n-1
    if J~=K
        res = res*sin(pi*(K-J)/(2*n))*sin(pi*(K+J)/(2*n));
    end
end 
res = 1/res;