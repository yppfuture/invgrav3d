function [m0]=conjgrad3D(m0,G,d,nX,nY,nZ)

%Compute residual iteratively. Computes in order:
% 1- the derivative terms WxtWx*m , WytWy*m, WztWz*m
% 2- GtWdtWdG*m
% 3- The smallest term WstWs*m
% 4- The right-hand-side G'*WdtWd*d
% r=(A*(x0)-RHS);

wd=std(d);
tic
LHS=CompLHS(m0,G,nX,nY,nZ);
RHS=G'*d;

r=(LHS-RHS);
p=-r;

% invert=zeros(20,60);
rold=r'*r;
rnorm=norm(r);
dnorm=norm(RHS);
count=1;
error=rnorm/dnorm;
misfit(count)=norm(r)^2;

toc
%     figure(5)
%     plot (count,misfit(count),'*');
%     hold on
while (count)<30%error>=10e-4
    tic
    Ap=CompLHS(p,G,nX,nY,nZ);
    alpha=rold./(p'*Ap);
   
    
     m0=m0+alpha.*p;

    r=r+alpha*Ap;       %+(0*drxz').^2);
    rnew=r'*r;
    p=-r+rnew/rold.*p;
    rold=rnew;

   
%     invert(:)=x0(:);
    
    rnorm=norm(r);
    error=rnorm/dnorm;
    count=count+1;
    misfit(count)=norm(r)^2;
    
%     figure(5)
%     plot (count,misfit(count),'*');
%     hold on
    
    
    toc
end

