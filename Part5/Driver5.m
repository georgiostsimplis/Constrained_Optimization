%=============Driver file for 5.5==================%
%==================="fmincon"======================%
%==================Initial data====================%
clc; clear;
LB=[-5;0]; UB=[-1;5];
x0=[-5;-5];
A=[]; B=[]; Aeq=[]; Beq=[];
%=================Call function===================%
options = optimoptions( 'fmincon','Display','iter',...
    'SpecifyObjectiveGradient',true,...
    'SpecifyConstraintGradient',true,'FunctionTolerance',10e-7);

[x,fval,exitflag,output]= fmincon(@(x) fun2(x),x0,A,B,Aeq,Beq,LB,UB,...
    @(x) con12(x),options)
%===========Specify obj function & gradient=======%
function [f,df] = fun2(x)
x1=x(1,1);
x2=x(2,1);

tmp1=x1^2+x2-11;
tmp2=x1+x2^2-7;

f=tmp1^2+tmp2^2;     %Gradient%
if nargout > 1
    df = zeros(2,1);
    df(1,1) = 4*tmp1*x(1) + 2*tmp2;
    df(2,1) = 2*tmp1 + 4*tmp2*x(2);
end
end
%============Specify constraint function=========%
function [c,ceq,dc,dceq]=con12(x)
x1=x(1,1);
x2=x(2,1);
tmp=x1+2;
ceq=tmp^2-x2;
c=[];
if nargout > 2
    dc =[];
    dceq = zeros(2,1);
    dceq(1,1) = 2*tmp;
    dceq(2,1) = -1.0; 
end
end

