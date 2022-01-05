function [x,lambda]=EqualityQPSolverLUsparse(H,g,A,b)
n=length(g);
m=length(b);
z1=zeros(m);
K=sparse([H -A;-A' z1]);
if det(K)~=0
    d=-[g;b];
    [L,U,p]=lu(K,'vector');
    z=U\(L\d(p));
    x=z(1:n);
    lambda=z(n+1:end);
else
    disp('There is not a unique minimizer')
end
end
