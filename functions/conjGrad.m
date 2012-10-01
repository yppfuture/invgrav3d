function m0 = conjGrad(m0, G, d, nx ,ny ,nz)

%Compute residual iteratively. Computes in order:
% 1- the derivative terms WxtWx*m , WytWy*m, WztWz*m
% 2- GtWdtWdG*m
% 3- The smallest term WstWs*m
% 4- The right-hand-side G'*WdtWd*d
% r=(A*(x0)-RHS);

%wd = std(d);

LHS = getLHS(m0, G, nx, ny, nz);
RHS = G' * d;

r = (LHS - RHS);
p = -r;

rold = r' * r;
%rnorm = norm(r);
%dnorm = norm(RHS);
count = 1;
%error = rnorm / dnorm;
%misfit(count) = norm(r)^2;

while (count) < 10%error>=10e-4

    Ap = getLHS(p, G, nx, ny, nz);
    alpha =rold ./ (p' * Ap);


    m0 = m0 + alpha .* p;

    r = r + alpha * Ap;
    rnew = r' * r;
    p = -r + rnew / rold .* p;
    rold = rnew;

    %rnorm = norm(r);
    %error = rnorm / dnorm;
    count = count + 1;
    %misfit(count) = norm(r)^2;


end
