%=============Line search method============%
function [xk,initer]=linesearch(x,y,z,p)
%=============Initialize data===============%
alpha=1; initer=0;
flag=0;
[f,df,~]=fun1(x);
[ceq,dceq,~]=con1(x);
[c,dc]=con2(x);
f0=f+y'*abs(ceq)+z'*abs(min(zeros(length(z),1),c));
f01=df'*p-y'*abs(ceq)-z'*abs(min(zeros(length(z),1),c));
while flag==0
%===Compute new point & functions' value===% 
    xk=x+alpha*p;
    [f,df,~]=fun1(xk);
    [ceq,dceq,~]=con1(xk);
    [c,dc]=con2(xk);
    fa=f+y'*abs(ceq)+z'*abs(min(zeros(length(z),1),c));
%===========Check armijo condition=========%
    if fa<=f0+0.1*f01*alpha
        flag=1;
    else
%============compute new alpha=============%
        a1=(fa-(f0+f01*alpha))/(alpha*alpha);
        amin=-f01/(2*a1);
        alpha=min(0.9*alpha,max(amin,0.1*alpha));
        initer=initer+1;
    end
end
end
        