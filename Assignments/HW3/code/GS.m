function [Q, R] = GS(A)
    %classical gram schmidt factorization on matrix A    
    [nRows, nCols] = size(A);
    Q = zeros(nRows, nCols);
    Q(:,1) = A(:,1)/norm(A(:,1));
    for k=2:nCols
       Q(:,k) = A(:,k);
       for n=1:k-1
           Q(:,k) = Q(:,k) - ((A(:,k)'*Q(:,n))/(Q(:,n)'*Q(:,n)))*Q(:,n);
       end
       Q(:,k) = Q(:,k)/norm(Q(:,k));
    end
    R = transpose(Q)*A;
end