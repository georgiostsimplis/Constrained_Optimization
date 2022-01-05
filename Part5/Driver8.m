%===========SQP Trust region algorithm========%
%===============Initialize data===============%
clc; clear;
x=[-3;-3]; y=1;
z1=ones(2,1); z2=ones(2,1);
k=0; kmax=50;
[f,df,d2f]=fun1(x); %functions at starting point
[ceq,dceq,d2ceq]=con1(x);
[c,dc]=con2(x);
lb=[-5;0];
ub=[-1;5];
b1=[0;5;5;-1];

H=d2f; %Hessian
tol=10e-7;
options =  optimoptions('quadprog','Display','off');
dk=inf; %Starting trust region
converged=false;
initer=0;
t1=tic;
%=====Enter loop to solve sub-problems=====%
while ~converged &&(k<kmax)
%=========Combine LB-UB constraints========%  
    lb1=lb-x(:,end);
    ub1=ub-x(:,end);
    lbb=-dk*ones(2,1);
    ubb=dk*ones(2,1);
    lb1=max(lb1,lbb);
    ub1=min(ub1,ubb);
    [p,~,~,~,lambda] = quadprog(H,df,[],[],dceq',-ceq,lb1,ub1,[],options);
    y=-lambda.eqlin;
    z1=lambda.lower;
    z2=lambda.upper;
   
    [f1,df1,d2f1]=fun1(x(:,end)+p);
%===========Conditions for ro============%
    if p~=0
        ro=(f1-f)/(0.5*p'*H*p+df'*p);
    else
        ro=1;
    end
%===========Conditions for gamma============%    
    if ro<0.25
        gamma=0.25;
    elseif (ro<=0.75) && (ro>=0.25)
        gamma=1;
    else
        gamma=2;
    end
%================Accept step================%
    if ro>0
        k=k+1;
        dk=gamma*dk;
        x(:,end+1)=x(:,end)+p;
        gL=df-dceq*y-dc(:,1)*z1(1)-dc(:,2)*z1(2)-dc(:,3)*z2(1)-dc(:,4)*z2(2);
        [f,df,d2f]=fun1(x(:,end));
        [ceq,dceq,d2ceq]=con1(x(:,end));
        [c,dc]=con2(x(:,end));
        gL1=df-dceq*y-dc(:,1)*z1(1)-dc(:,2)*z1(2)-dc(:,3)*z2(1)-dc(:,4)*z2(2);
        q=gL1-gL;
        H= Damped_BFGS(H,q,p);
        converged=(norm(gL1,inf)<tol) &&(norm(p,inf)<tol);
    else
%================Reject step================%
        dk=norm(p,inf)*gamma;
        initer=initer+1;
    end
end
time=toc(t1)
k
x
initer
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
plot(k1,x(1,:),'Linewidth',1.5)
plot(k1,x(2,:),'Linewidth',1.5)
hold off