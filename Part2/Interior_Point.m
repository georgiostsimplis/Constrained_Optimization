%         Convex Quadratic Programming
%     Primal-Dual Interior Point Algorithm
%===============================================
%             min (1/2)*x'*H*x + g'*x
%        s.t        A'*x = b
%                   C'*x >= d
%===============================================

function[x,fval,y,z,iter,t] = Interior_Point(H,g,A,b,C,d,x0)
tStart=tic;
n = length(H); [~,mc] = size(C); [~,n1] = size(A);
x = x0; y = ones(n1,1);z = ones(mc,1);s = ones(mc,1);
itermax = 150;
iter = 0;

tol = 10e-9;
%compute the residuals
rL = H*x+g-A*y-C*z;
rA = b-A'*x;
rC = d-C'*x+s;
S = diag(s);Z=diag(z);e=ones(mc,1);
rsz = S*Z*e;
mu = (z'*s)/mc;
 %Stopping Criteria   
Converged = (norm(rL, inf) < tol ) && (norm(rA , inf) < tol) &&...
    (abs(mu)<tol);
%Start the loop
while ~Converged && (iter<itermax)
    iter = iter+1;
%First System to affine direction
    H1=H+C*(inv(S)*Z)*C';
    KKT=[H1 -A;-A' zeros(n1)];
    rL1=rL-C*(inv(S)*Z)*(rC-inv(Z)*rsz);
    RH=-[rL1;rA];
    %LDL factorization
    [L,D,p] = ldl(KKT,'lower','vector');
    sol(p) = L'\(D\(L\RH(p)));
    %Affine direction
    dxaff = sol(1:n)'; dyaff = sol(n+1:n+n1)'; 
    dzaff=-(inv(S)*Z)*C'*dxaff+(inv(S)*Z)*(rC-inv(Z)*rsz);
    dsaff=-inv(Z)*rsz-inv(Z)*S*dzaff;
    idx = find(dzaff < 0.0);
    alpha = min([1.0; -z(idx,1)./dzaff(idx,1)]);
    
    idx = find(dsaff < 0.0);
    beta = min([1.0; -s(idx,1)./dsaff(idx,1)]);
    alpha_aff = min(alpha,beta);
   %Duality gap 
    muaff = ((z+alpha_aff*dzaff)'*(s+alpha_aff*dsaff))/mc;
    
    sigma = (muaff/mu)^3;
    %Affine-Centering-Correction Direction
    rsz1=rsz+diag(dsaff)*diag(dzaff)*e-sigma*mu*e;
    rL1=rL-C*(inv(S)*Z)*(rC-inv(Z)*rsz1);
    RH=-[rL1;rA];
    [L,D,p] = ldl(KKT,'lower','vector');
    sol(p) = L'\(D\(L\RH(p)));
% Solution to obtain direction toward optimal
    dx = sol(1:n)'; dy = sol(n+1:n+n1)'; 
    dz=-(inv(S)*Z)*C'*dx+(inv(S)*Z)*(rC-inv(Z)*rsz1);
    ds=-inv(Z)*rsz1-inv(Z)*S*dz;
    
    
    
    idx = find(dz < 0.0);
    alpha = min([1.0; -z(idx,1)./dz(idx,1)]);
    
    idx = find(ds < 0.0);
    beta = min([1.0; -s(idx,1)./ds(idx,1)]);
    alpha = min(alpha,beta); 
  %Move to the new point  
    x=x+0.995*alpha*dx;
    y=y+0.995*alpha*dy;
    z=z+0.995*alpha*dz;
    s=s+0.995*alpha*ds;
    
    rL = H*x+g-A*y-C*z;
    rA = b-A'*x;
    rC = d-C'*x+s;
    S = diag(s);Z = diag(z);e = ones(mc,1);
    rsz = S*Z*e;
    mu = (z'*s)/mc;
    
    Converged = (norm(rL, inf) < tol ) && (norm(rA , inf) < tol) &&...
        (abs(mu)<tol);
end
x = x(:,end);
fval = 0.5*x'*H*x + g'*x;
t = toc(tStart);
end