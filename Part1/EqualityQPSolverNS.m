function [x,lambda]=EqualityQPSolverNS(H,g,A,b)
n = length(g);
[Q,Rbar] = qr(A);
m1 = size(Rbar,2);
Q1 = Q(:,1:m1);
Q2 = Q(:,m1+1:n);
R = Rbar(1:m1,1:m1);
xy = R'\b;
xz = (Q2'*H*Q2)\(-Q2'*(H*Q1*xy+g));
x = Q1*xy+Q2*xz;
lambda = R\(Q1'*(H*x+g));
end