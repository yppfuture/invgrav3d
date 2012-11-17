function [m0,error] = conjGrad(m0, G, d, Wr, WctWc, wd, nx ,ny ,nz, dx, dy, dz, beta, alphaS)


%Compute residual iteratively. Computes in order:
% 1- the derivative terms WxtWx*m , WytWy*m, WztWz*m
% 2- GtWdtWdG*m
% 3- The smallest term WstWs*m
% 4- The right-hand-side G'*WdtWd*d
% r=(A*(x0)-RHS);

%wd = std(d);

tic
LHS = compLHS_v6(m0, G, Wr, WctWc, wd, nx, ny, nz, dx, dy, dz, beta, alphaS);
toc
RHS = G' * d;

r = (LHS - RHS);
p = -r;

rold = r' * r;
rnorm = norm(r);
dnorm = norm(RHS);
count = 1;
error(count) = rnorm / dnorm;
% figure(1)=plot(count,error);
% hold on
% misfit(count) = norm(r)^2;

while (count) < 5%error>=10e-4
tic

    Ap = compLHS_v6(p, G, Wr, WctWc, wd, nx, ny, nz, dx, dy, dz, beta, alphaS);

    alpha =rold ./ (p' * Ap);


    m0 = m0 + alpha .* p;

    r = r + alpha * Ap;
    rnew = r' * r;
    p = -r + rnew / rold .* p;
    rold = rnew;
    
    count = count + 1;
    rnorm = norm(r);
    error(count) = rnorm / dnorm;
    
%     figure(1)=plot(count,error);
%     hold on
    %misfit(count) = norm(r)^2;
toc

end
