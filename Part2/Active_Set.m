%             Active Set Algorithm
%             Quadratic Programming
%==================================================
% min (1/2)*x'*G*x + g'*x
%  s.t A'*x = b
%      C'*x >=d
%==================================================
% k = iterations
% t = cputime
function [x,lambda,mu,fval,t,k] = Active_Set(G,g,A,b,C,d,x0)
tStart = tic;
workingSet = find(C'*x0-d == 0);
nonworkingSet = find(C'*x0-d ~= 0);
flag = 0;
k = 0;
x = x0;
mu=zeros(length(C),1);
maxiter=100;
[~,n2]=size(A);


while flag == 0 && k<maxiter
    k = k+1;
    n1=length(workingSet);
    
   
    KKT = [G -A -C(:,workingSet);
           -A' zeros(n2) zeros(n2,n1);
           -C(:,workingSet)' zeros(n1,n2) zeros(n1)];
    RHS = -[G*x(:,k)+g; b-A'*x(:,k); d(workingSet)-C(:,workingSet)'*x(:,k)];
    
    sol = KKT\RHS;
    p = sol(1:length(G));
    lambda=sol(length(G)+1 : length(G)+n2);
    mu(workingSet) = sol(length(G)+n2+1 : end);
    
    if norm(p)<10e-9
        if all(mu(workingSet)>=0)
            flag = 1;
            k =  k-1;
            mu(nonworkingSet)=0;
        else
            [minimum,j1] = min(mu(workingSet));
            j=find(mu==minimum);
            nonworkingSet = sort([nonworkingSet;j]);
            workingSet(j1) = [];
            x(:,k+1) = x(:,k);
        end
    else
        temp1 = C(:,nonworkingSet)'*p;
        index = find(temp1<0);
        i = nonworkingSet(index);
        temp = (d(i)-C(:,i)'*x(:,k))./(C(:,i)'*p);
        [alpha,j] = min(temp);
        j = i(j);
        if alpha<1
            x(:,k+1) = x(:,k)+alpha*p;
            workingSet = sort([workingSet;j]);
            j = find(nonworkingSet==j);
            nonworkingSet(j) = [];
        else
            x(:,k+1) = x(:,k)+p;
        end
    end
end
fval=0.5*x(:,end)'*G*x(:,end) + g'*x(:,end);
x=x(:,end);
t=toc(tStart);
mu;
end