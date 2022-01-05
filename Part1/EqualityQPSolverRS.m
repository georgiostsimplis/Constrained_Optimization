function [x,lambda]=EqualityQPSolverRS(H,g,A,b)
if det(H)~=0
    L = chol(H,'lower');
    v = L'\(L\g);
    HA = A'*H^-1*A;
    lambda = HA\(b+A'*v);
    x = H\(A*lambda-g);
else
    disp('There is not a unique minimizer')
end
end