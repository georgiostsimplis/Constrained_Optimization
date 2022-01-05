%=====================================%
%         Driver file for 2.8         %
%=====================================%
clc; clear;
%=====Initialize data=================%
G = [2.3 0.93 0.62 0.74 -0.23;
     0.93 1.4 0.22 0.56 0.26;
     0.62 0.22 1.8 0.78 -0.27;
     0.74 0.56 0.78 3.4 -0.56;
     -0.23 0.26 -0.27 -0.56 2.6];
A = [15.1 1;
    12.5 1;
    14.7 1;
    9.02 1;
    17.68 1];
g=zeros(5,1); C=eye(5); b = [10; 1]; d=zeros(5,1); x0=0.2*ones(5,1);
%======Create matrix for results========%
matrix=cell(4,7);
matrix(2,1)={"Active-Set Algorithm"};
matrix(3,1)={"Interior-Point Algorithm"};
matrix(4,1)={"quadprog library"};
matrix(1,2)={"Optimal"};
matrix(1,3)={"lambda_Eq"};
matrix(1,4)={"lambda_Ineq"};
matrix(1,5)={"f-value"};
matrix(1,6)={"time"};
matrix(1,7)={"Iterations"};

%======Solution Active-Set & store data=====%
disp("==========Active-Set Algorithm============")
[x,lambda,mu,fval,t,k] = Active_Set(G,g,A,b,C,d,x0)
matrix(2,2)={x};
matrix(2,3)={lambda};
matrix(2,4)={mu};
matrix(2,5)={fval};
matrix(2,6)={t};
matrix(2,7)={k};
%======Solution Interior-Point & store data=====%
disp("========Interior-Point Algorithm=============")
H=G;
[x,fval,y,z,iter,t] = Interior_Point(H,g,A,b,C,d,x0)
matrix(3,2)={x};
matrix(3,3)={y};
matrix(3,4)={z};
matrix(3,5)={fval};
matrix(3,6)={t};
matrix(3,7)={iter};
%======Solution quadprog & store data=====%
disp("===========quadprog library=============")
f=g; Aeq=A';  beq=b;  lb=zeros(5,1); ub=ones(5,1);
t1=tic;
[x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,x0)
time=toc(t1)
matrix(4,2)={x};
matrix(4,3)={lambda.eqlin};
matrix(4,4)={'-'};
matrix(4,5)={fval};
matrix(4,6)={time};
matrix(4,7)={output.iterations};
disp("For summary solution click on matrix at the Workspace!!!!")

matrix