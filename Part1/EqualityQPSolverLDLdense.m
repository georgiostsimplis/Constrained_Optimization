function [x,lambda]=EqualityQPSolverLDLdense(H,g,A,b)
n=length(g);
m=length(b);
z1=zeros(m);
K=[H -A;-A' z1];
d=-[g;b];
if det(K)~=0
    z=zeros(m+n,1);
    [L,D,p]=ldl(K,'lower','vector');
    z(p)=L'\(D\(L\d(p)));
    x=z(1:n);
    lambda=z(n+1:end);
else
    disp('There is not a unique minimizer')
end
end