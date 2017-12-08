% adaptive quadrature using Simpson's rule
%
% input: f - integrand, scalar function that takes one input
% a,b - integration is over interval [a,b]
% ep - error tolerance
% count - number of times the integrand is evaluated
% this should not be passed in by the user, it
% is only used in recursive calls
%
% output: q - the approximation to the integral
% new_count - number of times the integrand was evaluated
%
% note: this function requires that the quadrature rule S be
% defined in the same file
%
function [q,new_count]=adapt(f,a,b,ep,count)
    % if count was not passed in, initialize it to zero
    %
    if( nargin < 5 )
        count = 0;
    end
    new_count = count + 9;
    sab = S(f,a,b);
    c = 0.5*(a+b);
    sac = S(f,a,c);
    scb = S(f,c,b);
    factor = 15;
    if( abs(sac+scb-sab) < factor*ep )
        q = sac + scb;
    else
        [q1,new_count]=adapt(f,a,c,ep/2,new_count);
        [q2,new_count]=adapt(f,c,b,ep/2,new_count);
        q =q1+q2;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simposon's Rule
% input: f - integrand, scalar function that takes one input
% a,b - integration is over [a,b]
%
% output: q - simpson's rule approximation to integral of f over [a,b]
%
function q=S(f,a,b)
    q = (b-a)/6*( feval(f,a) + 4*feval(f,(a+b)/2) + feval(f,b));
end
