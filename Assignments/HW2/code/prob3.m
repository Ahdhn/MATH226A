clc;
clear;
close all;
disp('Math 226A - HW2 - Problem #3 (b):');

m=2;
i=1;
while true
    A = randn(m);
    [L,U] = lu(A);
    my_rho =  max(max(abs(U))) / max(max(abs(A)));
    rhos(i,1)=m;
    rhos(i,2) = my_rho;    
    m=m+10;
    i = i+1;
    if m >= 5e3
        break
    end
end

P=polyfit(rhos(:,1),rhos(:,2),3);

scatter(rhos(:,1),rhos(:,2));
%set(gca,'xscale','log');
%set(gca,'yscale','log');
hold on
pol = polyval(P,rhos(:,1));
%plot(rhos(:,1),pol);
ylabel('Growth Factor','FontSize',15);
xlabel('Matrix Size','FontSize',15);
%set(gca,'xscale','log');
%set(gca,'yscale','log');

polyval(P,1e6)