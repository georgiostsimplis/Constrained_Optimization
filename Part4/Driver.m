%======Driver File 4.4, 4.6, 4.7==============
%       min f(x)=2*x1-5*x2
%        s.t    x1+x2=210
%               100<=x1<=200
%               50<=x2<=150
%=============================================
clc; clear;
A=[1 1 0 1 0;
    1 0 1 0 1
    0 -1 0 0 0;
    0 0 -1 0 0;
    0 0 0 1 0;
    0 0 0 0 1]';
b=[210 100 50 200 150]'; g=[2;-5;0;0;0;0];
disp("========Simplex Algorithm=========")
t1=tic;
[x,mu,lambda,optimal,iter]=Simplex(g,A,b)
time1=toc(t1)
disp("========Interior-Point Algorithm=========")
x=ones(6,1);
t2=tic;
[x,info,mu,lambda,iter] = LPippd(g,A,b,x)
time2=toc(t2)
disp("============linprog===========")
f=g(1:2); Aeq=[1 1]; beq=210; lb=[100;50]; ub=[200;150];
t3=tic;
[x,fval,exitflag,output,lambda] = linprog(f,[],[],Aeq,beq,lb,ub)
time3=toc(t3)
