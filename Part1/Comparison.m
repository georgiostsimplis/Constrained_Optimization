%==============Driver file for 1.4====================%
clc;clear;
%====Create symmetric positive definite matrix H======%
H=randi([1,3],6);
H=H*H';
%========Create random column vectors g, A, b=========%
g=randi([1,9],6,1);
A=randi([1,9],6,1);
b=randi([1,9]);
k=zeros(5,1);
%=========for loop to solve in sequence QP's==========%
for i=1:5
    %=======Check-Create A with full column rank======% 
   [U,S,V]=svd(A);
   ra=i;
   A=U(:, 1: ra)* S(1: ra, 1: ra)* V(:, 1: ra)';
   k(i,1)=i;
   %==========Solve the QP's and store time==========%
   tic;
   [x,lambda]=EqualityQPSolver(H,g,A,b,'LUdense');
   t(i,1)=toc;
   tic;
   [x,lambda]=EqualityQPSolver(H,g,A,b,'LUsparse');
   t(i,2)=toc;
   tic;
   [x,lambda]=EqualityQPSolver(H,g,A,b,'LDLdense');
   t(i,3)=toc;
   tic;
   [x,lambda]=EqualityQPSolver(H,g,A,b,'LDLsparse');
   t(i,4)=toc;
   tic;
   [x,lambda]=EqualityQPSolver(H,g,A,b,'NullSpace');
   t(i,5)=toc;
   tic;
   [x,lambda]=EqualityQPSolver(H,g,A,b,'RangeSpace');
   t(i,6)=toc;
   %===========Add another column constraint=========%
   A1=randi([1,9],6,1);
   b1=randi([1,9]);
   A=[A A1];
   b=[b;b1];
end
%===============Plot the results=====================%
plot(k,t(:,1),'r',k,t(:,2),'m',k,t(:,3),'y',k,t(:,4),'c',k,t(:,5),...
    'b',k,t(:,6),'g','Linewidth',2)
lgd=legend('LUdense','LUsparse','LDLdense','LDLsparse','NullSpace',...
    'RangeSpace');
lgd.NumColumns = 2;
title("Solvers' behavior in respect to number of Constraints");
xticks([1,2,3,4,5])
xlabel('Number of constraints');
ylabel('Time');