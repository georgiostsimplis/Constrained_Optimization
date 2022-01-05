%========Driver file for 3.5=======%
clc;clear;
H = [2.3 0.93 0.62 0.74 -0.23 0;
     0.93 1.4 0.22 0.56 0.26 0;
     0.62 0.22 1.8 0.78 -0.27 0;
     0.74 0.56 0.78 3.4 -0.56 0;
     -0.23 0.26 -0.27 -0.56 2.6 0;
     0 0 0 0 0 0];
Aeq = [15.1 1;
    12.5 1;
    14.7 1;
    9.02 1;
    17.68 1;
    2 1]';
f=zeros(6,1); x0=[0.1 0.1 0.2 0.2 0.2 0.2]';
lb=zeros(6,1); ub=ones(6,1);

risk=zeros(10,2);

i=1;
options = optimoptions('quadprog','Display','off');
for R=2:0.2:17.6
    beq=[R;1];
    [x,fval] = quadprog(H,f,[],[],Aeq,beq,lb,ub,x0,options);
    risk(i,1)=R;
    %==== Compute Variance=========%
    risk(i,2)=2*fval;
    %=========== Store data ===============%
    sol1(1,i)=x(1);sol1(2,i)=R;
    sol2(1,i)=x(2);sol2(2,i)=R;
    sol3(1,i)=x(3);sol3(2,i)=R;
    sol4(1,i)=x(4);sol4(2,i)=R;
    sol5(1,i)=x(5);sol5(2,i)=R;
    sol6(1,i)=x(6);sol6(2,i)=R;
    i=i+1;
end

%====Efficient frontier plot with securities===%
figure(1)
hold on
plot(risk(:,2),risk(:,1),'LineWidth',2)
legendInfo{1} = ['Efficient frontier'];

s=diag(H);
for j=1:6
   plot(s(j),Aeq(1,j),'X','LineWidth',3)
   legendInfo{j+1} = ['Security ' num2str(j)];
end
legend(legendInfo)
title("Efficient Frontier")
xlabel("Variance")
ylabel("Expected Return")
hold off
%============ optimal portfolio=============%
figure(2)
hold on 
plot(sol1(2,:),sol1(1,:),'-','LineWidth',2,'DisplayName','security 1')
plot(sol2(2,:),sol2(1,:),'-','LineWidth',2,'DisplayName','security 2')
plot(sol3(2,:),sol3(1,:),'-','LineWidth',2,'DisplayName','security 3')
plot(sol4(2,:),sol4(1,:),'-','LineWidth',2,'DisplayName','security 4')
plot(sol5(2,:),sol5(1,:),'-','LineWidth',2,'DisplayName','security 5')
plot(sol6(2,:),sol6(1,:),'-','LineWidth',2,'DisplayName','free risk sec.')
plot(sol1(2,41),sol1(1,41),'X','LineWidth',2)
plot(sol2(2,41),sol2(1,41),'X','LineWidth',2)
plot(sol3(2,41),sol3(1,41),'X','LineWidth',2)
plot(sol4(2,41),sol4(1,41),'X','LineWidth',2)
plot(sol5(2,41),sol5(1,41),'X','LineWidth',2)
plot(sol6(2,41),sol6(1,41),'X','LineWidth',2)

legend('security 1','security 2','security 3','security 4','security 5','free risk sec.')
ylabel("Portfolio's weighting  factors ")
xlabel("Return")
hold off
%========optimal portfolio results for R=10=========%
H = [2.3 0.93 0.62 0.74 -0.23 0;
     0.93 1.4 0.22 0.56 0.26 0;
     0.62 0.22 1.8 0.78 -0.27 0;
     0.74 0.56 0.78 3.4 -0.56 0;
     -0.23 0.26 -0.27 -0.56 2.6 0;
     0 0 0 0 0 0];
Aeq = [15.1 1;
    12.5 1;
    14.7 1;
    9.02 1;
    17.68 1;
    2 1]';
f=zeros(6,1); x0=[0.1 0.1 0.2 0.2 0.2 0.2]';
lb=zeros(6,1); ub=ones(6,1);
R=10; beq=[R;1];
disp("==== Optimal Portfolio with risk free security ====")
[x,fval] = quadprog(H,f,[],[],Aeq,beq,lb,ub,x0,options)

