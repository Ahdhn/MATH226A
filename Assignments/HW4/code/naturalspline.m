% Natural spline coefficients
%
% Input: x - n by 1 vector of the nodes
% y - n by 1 vector with the corresponding function values
% Output: P - n-1 by 4 matrix
% the elements of the ith row of P give the coefficients
% of the cubic between [x(i),x(i+1)] as
% Si(x) = P(i,1) + P(i,2)*(x-x(i))
% + P(i,3)*(x-x(i))^2 + P(i,4)*(x-x(i))^3
%
function P=naturalspline(x,y,not_a_knot)
    % initialize P
    %    
    n = length(x);
    P = zeros(n-1,4);
    % compute the distances between the points
    %
    h = x(2:n) - x(1:n-1);
    % set up the tridiagonal linear system to solve
    %
    d0 = 2*(h(1:n-2)+h(2:n-1)); % diagonal
    d1 = h(2:n-2); % super and sub diagonal
    % form the matrix and rhs
    %
    A = diag(d0,0) + diag(d1,-1) + diag(d1,1);
    b = 3*(y(3:n)-y(2:n-1))./h(2:n-1) - 3*(y(2:n-1)-y(1:n-2))./h(1:n-2);
    
    if not_a_knot 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %New Conditions for not-a-knot
        A = [zeros(1,n-2) ; A ; zeros(1,n-2)];
        A = [zeros(n,1), A, zeros(n,1)];  
        %From the 1st boundary condition 
        A(1,1)=h(2);%c_{0}
        A(1,2)=-h(2)-h(1);%c_{1}
        A(1,3)=h(1);%c_{2}
        
        %From the original equation at j=1
        A(2,1)=h(1);%c_{0}
        A(2,2)=2*(h(1)+h(2));%c_{1}
        A(2,3)=h(2);%c_{2}
        %From the original equation at j=n-1
        A(n-1,n-2)=h(n-2); %
        A(n-1,n-1)=2*(h(n-2)+h(n-1));
        A(n-1,n)=h(n-1);
        
        %From the 2nd boundary condition 
        A(n,n-2)=h(n-1); %c_{n-2}
        A(n,n-1)=-h(n-2)-h(n-1);%c_{n-1}
        A(n,n)=h(n-2);%c_{n}
        
        %%% Update the right hand side 
        b =[0;b;0];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    % solve the linear system
    %
    z = A\b;
    % append the natural boundary conditions
    %
    if ~not_a_knot
        z = [0; z; 0];
    end
    % compute the coefficients
    %
    P(:,1) = y(1:n-1);
    P(:,2) = (y(2:n)-y(1:n-1))./h - h.*(z(2:n)+2*z(1:n-1))/3;
    P(:,3) = z(1:n-1);
    P(:,4) = (z(2:n)-z(1:n-1))./(3*h);
end