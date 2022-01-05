%   Revised Simplex Algorithm
%======================================
%  LP in Standard form 
%      min f(x)= g'*x
%        s.t   A*x=b
%======================================


function [x,mu,lambda,optimal,iter]=Simplex(g,A,b)
%=====Initialize Data=====================
[m,n] = size(A);
basic=1;
lambda = zeros(n,1);
iter = 0;
set=[1:n];
x=-1*ones(n,1);
flag = 0;
%slack=find(g==0)
%=====Choose Basic and Non Basic Sets=====
while ~all(x(basic)>=0)
    non_basic = randperm(n,n-m);
    basic = setdiff(set,non_basic);
    B = A(:,basic);
    N = A(:,non_basic);
    if det(B)~=0
        x(basic) = pinv(B)*b;
        x(non_basic)=0;
        mu = B'\g(basic);
        lambda(non_basic) = g(non_basic)-N'*mu;
        lambda(basic)=0;
    else
        x=-1*ones(n,1);
    end
end

%=====Simplex Loop=======
while flag==0
    iter = iter + 1;
    B = A(:,basic);
    N = A(:,non_basic);
    x(basic) = B\b;
    x(non_basic)=0;
    mu = B'\g(basic);
    lambda(non_basic) = g(non_basic)-N'*mu;
    lambda(basic)=0;
%=====Condition for break=====
    if lambda(non_basic)>=0
        flag = 1;
        optimal = g(basic)'*x(basic);
    else
%=====Conditions to change Sets====
        [lambda_min,s] = min(lambda(non_basic));
        s = find(lambda==lambda_min);
        h = B\A(:,s);
        index = find(h>0);
        [~,J] = min(x(basic(index))./h(index));
        J = index(J);
        J = basic(J);
        if isempty(J)
            flag = 1;
            disp("*** The LP is Unbounded ***")
            x=[]; mu=[]; lambda=[]; optimal=[];
        else
%======Change The Sets======
            non_basic(end+1) = J;
            basic(end+1) = s;
            non_basic = non_basic(non_basic ~= s);
            basic = basic( basic ~= J);
        end
    end
end
