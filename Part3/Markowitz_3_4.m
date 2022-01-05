%============Driver file for 3.3-3.4==========%
clc; clear;
H = [2.3 0.93 0.62 0.74 -0.23;
     0.93 1.4 0.22 0.56 0.26;
     0.62 0.22 1.8 0.78 -0.27;
     0.74 0.56 0.78 3.4 -0.56;
     -0.23 0.26 -0.27 -0.56 2.6];
Aeq = [15.1 1;
    12.5 1;
    14.7 1;
    9.02 1;
    17.68 1]';
f=zeros(5,1); x0=0.2*ones(5,1);
lb=zeros(5,1); ub=ones(5,1);

risk=zeros(10,2);

i=1;
%==============collect data to plot=============%
for R=9:0.2:17.6
    beq=[R;1];
    options = optimoptions('quadprog','Display','off');
    [x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,x0,options);
    if R==10
        disp(x)
        disp(2*fval)
    end
    risk(i,1)=R;
    risk(i,2)=2*fval;
    sol1(1,i)=x(1);sol1(2,i)=R;
    sol2(1,i)=x(2);sol2(2,i)=R;
    sol3(1,i)=x(3);sol3(2,i)=R;
    sol4(1,i)=x(4);sol4(2,i)=R;
    sol5(1,i)=x(5);sol5(2,i)=R;
    i=i+1;
end
figure(1)
plot(risk(:,2),risk(:,1),'-','LineWidth',2)
title("Efficient Frontier")
xlabel("Variance")
ylabel("Expected Return")

figure(2)
hold on 
plot(sol1(2,:),sol1(1,:),'-','LineWidth',2,'DisplayName','security 1')
plot(sol2(2,:),sol2(1,:),'-','LineWidth',2,'DisplayName','security 2')
plot(sol3(2,:),sol3(1,:),'-','LineWidth',2,'DisplayName','security 3')
plot(sol4(2,:),sol4(1,:),'-','LineWidth',2,'DisplayName','security 4')
plot(sol5(2,:),sol5(1,:),'-','LineWidth',2,'DisplayName','security 5')
legend('security 1','security 2','security 3','security 4','security 5')
ylabel("Portfolio's weighting  factors ")
xlabel("Return")
hold off

