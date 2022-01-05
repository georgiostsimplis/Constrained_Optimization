%=============Driver File for 1.4=============%
%==============Input matrices=================%
clc;clear;
H=[6 2 1;2 5 2;1 2 4];
g=[-8;-3;-3];
A=[1 0;0 1;1 1];
b=[3;0];
%==============Switch among solvers===========%
disp("=========LU-Dense=======")
[x,lambda]=EqualityQPSolver(H,g,A,b,'LUdense')
disp("=========LU-Sparce=======")
[x,lambda]=EqualityQPSolver(H,g,A,b,'LUsparse')
disp("=========LDL-Dense=======")
[x,lambda]=EqualityQPSolver(H,g,A,b,'LDLdense')
disp("=========LDL-Sparse=======")
[x,lambda]=EqualityQPSolver(H,g,A,b,'LDLsparse')
disp("=========NullSpace=======")
[x,lambda]=EqualityQPSolver(H,g,A,b,'NullSpace')
disp("=========RangeSpace=======")
[x,lambda]=EqualityQPSolver(H,g,A,b,'RangeSpace')
