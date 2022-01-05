function[f,df,d2f]=fun1(x)
x1=x(1,1);
x2=x(2,1);

tmp1=x1^2+x2-11;
tmp2=x1+x2^2-7;

f=tmp1^2+tmp2^2;
df=zeros(2,1);
df(1,1)=4*x1*tmp1+2*tmp2;
df(2,1)=2*tmp1+4*x2*tmp2;

d2f=zeros(2,2);
d2f(1,1)=4*tmp1+8*x1^2+2;
d2f(2,1)=4*(x1+x2);
d2f(1,2)=d2f(2,1);
d2f(2,2)=4*tmp2+8*x2^2+2;