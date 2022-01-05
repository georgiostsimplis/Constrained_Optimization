%==========SQP Algorithm with Damped BFGS=======%

%===============Initialize data================%
clc; clear;
x=[-5;-5]; y=1;
%Lagange multipliers for lower-upper bound constraints
z1=ones(2,1); z2=ones(2,1);
k=0;%Iterations
kmax=80;
lb=[-5;0]; ub=[-1;5];
%================call functions================%
[f,df,d2f]=fun1(x);
[ceq,dceq,d2ceq]=con1(x);
[c,dc]=con2(x);
options =  optimoptions('quadprog','Display','off');
H=eye(2); %Initial Hessian
tol=10e-7;
converged=false; t1=tic;
while ~converged && k<kmax
    k=k+1;
 %=========Set the lower and upper bounds=======%
    lb1=lb-x(:,end);
    ub1=ub-x(:,end);
 %====================Solve QP==================%   
    [p,~,~,~,lambda] = quadprog(H,df,[],[],dceq',-ceq,lb1,ub1,[],options);
 %===============Obtain multipliers=============%
    y=-lambda.eqlin;
    z1=lambda.lower;
    z2=lambda.upper;
 %===============Gradient of Lagrangian=========%   
    dL=df-dceq*y-dc(:,1)*z1(1)-dc(:,2)*z1(2)-dc(:,3)*z2(1)-dc(:,4)*z2(2);
    x(:,k+1)=x(:,k)+p; %take step
 %==========Function values at new point========%  
    [f,df,d2f]=fun1(x(:,end));
    [ceq,dceq,d2ceq]=con1(x(:,end));
    [c,dc]=con2(x(:,end));
    
    dL1=df-dceq*y-dc(:,1)*z1(1)-dc(:,2)*z1(2)-dc(:,3)*z2(1)-dc(:,4)*z2(2);
    q=dL1-dL;
%=======Call Damped BFGS to update Hessian======%
    H= Damped_BFGS(H,q,p);
%===========Set convergence criteria============%
    converged= ( norm(dL1,inf)<tol) && (norm(ceq,inf)<tol);
    
end
time=toc(t1)
x(:,end)
k
%%
%=============PLOT 1- point movements===========%
figure(1)
X1=-6:0.005:4;
X2=-5:0.005:5;
[X1,X2]=meshgrid(X1,X2);
f=(X1.^2+X2-11).^2+(X1+X2.^2-7).^2;
v = [0:2:10 10:10:100 100:40:1000];
contour(X1,X2,f,v,'LineWidth',1.5)
xlabel('x_{1}')
ylabel('x_{2}')
zlabel('f')

x2=-6:0.01:4;
y2=(x2+2).^2;
hold on 
plot(x2,y2,'k','LineWidth',3)
xlim([-6 4])
ylim([-5 5])

fill([X1(1) -5 -5 X1(1)],[X2(1) X2(1) X2(end) X2(end)],...
    'k','facealpha',0.2)
fill([-1 X1(end) X1(end) -1],[X2(1) X2(1) X2(end) X2(end)],...
    'k','facealpha',0.2)
fill([X1(1) X1(end) X1(end) X1(1)],[X2(1) X2(1) 0 0],...
    'k','facealpha',0.2)
fill([X1(1) X1(end) X1(end) X1(1)],[5 5 X2(end) X2(end)],...
    'k','facealpha',0.2)
hold off

hold on 
plot(x(1,:),x(2,:),'Linewidth',1.5)
hold off
%=============Convergence of the point===========%
k1=[0:k];
figure(2)
hold on
plot(k1,x(1,:))
plot(k1,x(2,:))
hold off