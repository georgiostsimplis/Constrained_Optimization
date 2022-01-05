%=============Driver File for 4.8===============%
%=================Insert Data===================%
clc
clear
g=-1*[15.10;12.50;14.70;9.02;17.68];
A=[1 1 1 1 1];
b=[1]; x=ones(5,1);
disp("========Interior-Point Algorithm=========")
[x,info,mu,lambda,iter] = LPippd(g,A,b,x)
disp("========Simplex Algorithm=========")
[x,mu,lambda,optimal,iter]=Simplex(g,A,b)

Aeq=[1 1 1 1 1];
beq=1;
lb=zeros(5,1);
ub=[];
f=g;
disp("========'linpgrog'=========")
[x,fval,exitflag,output,lambda] = linprog(f,[],[],Aeq,beq,lb,ub)
