function H= Damped_BFGS(H,q,p);

if p'*q>=0.2*p'*(H*p)
    theta=1;
else
    theta=(0.8*p'*(H*p))/(p'*(H*p)-p'*q);
end
ro=theta*q+(1-theta)*(H*p);

H=H-(((H*p)*(H*p)')/(p'*(H*p)))+(ro*ro')/(p'*ro);
end