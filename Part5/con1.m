function [ceq,dceq,d2ceq]=con1(x)
x1=x(1,1);
x2=x(2,1);

tmp=x1+2;

ceq=tmp^2-x2;

dceq=zeros(2,1);
dceq(1,1)=2*tmp;
dceq(2,1)=-1;

d2ceq=zeros(2,2);
d2ceq(1,1)=2;