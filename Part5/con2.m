function [c,dc]=con2(x)
x1=x(1,1);
x2=x(2,1);

c=zeros(4,1);
c(1,1)=x1+5;
c(2,1)=x2;
c(3,1)=-x1-1;
c(4,1)=-x2+5;

dc=[eye(2) -eye(2)];