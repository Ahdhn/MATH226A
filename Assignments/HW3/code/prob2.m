%solving a system of equations using Newton's method with the option to
%include line search. 
clc;
clear;
close all;
disp('Math 226A - HW3 - Problem #2:');

%%%%%%%%%%%%%%%%%%% Part a)
global x11s x12s j1s f11s f12s; 
syms x11s x12s j1s f11s f12s;
f11s = (5-2*x11s)/(2*x12s-3);%function in symbolic 
f12s = (5-2*x12s)/(2*x11s-3);%to get the jacobian 
j1s = [jacobian(f11s,[x11s, x12s]);jacobian(f12s,[x11s, x12s])];
j1s = inv(j1s);%the inverse

disp('Results for part a) x0=(1 1), no Line Search');
[normF11, X11] = Newton(@Func1, @myJac1, [1 1]', 10 ,false);  %initial guess (1,1), no LS

disp('Results for part a) x0=(1 1), with Line Search');
[normF11LS, X11LS] = Newton(@Func1, @myJac1, [1 1]', 10 ,true);   %initial guess (1,1), LS

disp('Results for part a) x0=(3 3), no Line Search');
[normF12, X12] = Newton(@Func1, @myJac1, [3 3]', 10 ,false);  %initial guess (3,3), no LS

disp('Results for part a) x0=(3 3), with Line Search');
[normF12LS, X12LS] = Newton(@Func1, @myJac1, [3 3]', 10 ,true);   %initial guess (3,3), LS

disp('Results for part a) x0=(10 10), no Line Search');
[normF13, X13] = Newton(@Func1, @myJac1, [10 10]', 10 ,false);%initial guess (10,10), no LS

disp('Results for part a) x0=(10 10), with Line Search');
[normF13LS, X13LS] = Newton(@Func1, @myJac1, [10 10]', 10 ,true); %initial guess (10,10), LS

plotNorm(1,0,@Func1);
plot(X11(:,2),X11(:,1), 'c*');
plot(X11LS(:,2),X11LS(:,1), 'yo');

plotNorm(2,0,@Func1);
plot(X12(:,2),X12(:,1), 'c*');
plot(X12LS(:,2),X12LS(:,1), 'yo');

plotNorm(3,0,@Func1);
plot(X13(:,2),X13(:,1), 'c*');
plot(X13LS(:,2),X13LS(:,1), 'yo');
%%%%%%%%%%%%%%%%%%% Part b)
global x21s j2s f2s; 
syms x21s j2s f2s; 
f2s = (x21s^2)/(1+x21s^2);%function in symbolic 
j2s = jacobian(f2s,x21s);
j2s= inv(j2s);

disp('Results for part b) x0=1, no Line Search');
[normF21, X21] = Newton(@Func2, @myJac2, 1, 10, false); %initial guess (1), no LS

disp('Results for part b) x0=1, with Line Search');
[normF21LS, X21LS] = Newton(@Func2, @myJac2, 1, 10, true); %initial guess (1),  LS

disp('Results for part b) x0=0, no Line Search');
[normF22, X22] = Newton(@Func2, @myJac2, 10, 10, false);%initial guess (0), no LS

disp('Results for part b) x0=0, with Line Search');
[normF22LS, X22LS] = Newton(@Func2, @myJac2, 10, 10, true); %initial guess (0), LS

%%%%%%%%%%%%%%%%%%% Part c)
global x31s x32s j3s f31s f32s; 
syms x31s x32s j3s f31s f32s; 
f31s = x31s*x32s + x31s^3 +4;
f32s = x31s*x32s^2 + x32s +6;
j3s = [jacobian(f31s,[x31s, x32s]);jacobian(f32s,[x31s, x32s])];
j3s = inv(j3s);
 
disp('Results for part c) x0=(-1, -1), no Line Search');
[normF31, X31] = Newton(@Func3, @myJac3, [-1 -1]', 10 ,false);%initial guess (-1,-1), no LS

disp('Results for part c) x0=(-1, -1), with Line Search');
[normF31LS, X31LS] = Newton(@Func3, @myJac3, [-1 -1]', 10 ,true); %initial guess (-1,-1), LS

disp('Results for part c) x0=(1, 1), no Line Search');
[normF32, X32] = Newton(@Func3, @myJac3, [1 1]', 20 ,false);  %initial guess (1,1), no LS

disp('Results for part c) x0=(1, 1), with Line Search');
[normF32LS, X32LS] = Newton(@Func3, @myJac3, [1 1]', 20 ,true);   %initial guess (1,1), LS

plotNorm(4,0,@Func3);
plot(X31(:,2),X31(:,1), 'c-*');
plot(X31LS(:,2),X31LS(:,1), 'r:o');

plotNorm(5,0,@Func3);
plot(X32(:,2),X32(:,1), 'c-*');
plot(X32LS(:,2),X32LS(:,1), 'r:o');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%Plot Contours %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotNorm(fignum, X, func)
    v=[.25,.5:2:40];
    xr=-5:0.2:5;
    n=length(xr);
    z=zeros(n,n);
    for i=1:n
        for j=1:n        
            z(i,j)=norm(func([xr(i),xr(j)]));
        end
    end
    figure(fignum);
    contourf(xr,xr,z,v);
    ylim([-5 5]);
    xlim([-5 5]);
    hold 
    
   % plot(X(:,1),X(:,2),'-*');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%Jacobians & Functions%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Part a)
function j = myJac1(x)
    global x11s x12s j1s
    x11s = x(1);
    x12s = x(2);       
    j = double(subs(j1s));
end
function val = Func1(x)
    global x11s x12s f11s f12s
    %evluate the function on x
    %x is a vector
    x11s = x(1);
    x12s = x(2);
    val(1) = subs(f11s);
    val(2) = subs(f12s);
    val = double(val)';
end
%%%%%%%%%%%%%%%%%%% Part b)
function j = myJac2(x)
    global x21s j2s
    x21s = x;    
    j = double(subs(j2s));
end
function val = Func2(x)
    global x21s f2s   
    x21s = x;    
    val = double(subs(f2s));
end
%%%%%%%%%%%%%%%%%%% Part c)
function j = myJac3(x)
    global x31s x32s j3s
    x31s = x(1);
    x32s = x(2);
    j = subs(j3s);
    j = double(j);
end
function val = Func3(x)
    global x31s x32s f31s f32s
    %evluate the function on x
    %x is a vector
    x31s = x(1);
    x32s = x(2);
    val(1) = subs(f31s);
    val(2) = subs(f32s);
    val = double(val');
end


