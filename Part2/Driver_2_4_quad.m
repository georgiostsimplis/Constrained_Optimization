%============Driver File for 2.4 2.6 2.7============%
%====================Input matrices=================%
clc;clear;
A=[1;1]; b=0;
C=[1 0 -1 0;
  0 1 0 -1];
d=[-3 -3 -4 -4]';
H=[2,0;
   0,2];
g=[-2;-4]; x0=[0;0];
disp("===========Interior-Point Algorithm==========")
[x,fval,y,z,iter,t] = Interior_Point(H,g,A,b,C,d,x0)
G=H;
disp("===========Active-Set Algorithm==========")
[x,lambda,mu,fval,t,k] = Active_Set(G,g,A,b,C,d,x0)

Aeq=[1 1]; beq=0;
H=[2 0;0 2]; f=[-2;-4];
lb=[-3;-3]; ub=[4;4];
x0=[0;0];
t=tic;
disp("=============quadprog==============")
[x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,x0)
time=toc(t)