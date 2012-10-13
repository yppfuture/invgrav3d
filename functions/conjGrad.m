function m0 = conjGrad(m0, G, d, Wr, wd, nx ,ny ,nz, beta)

%Compute residual iteratively. Computes in order:
% 1- the derivative terms WxtWx*m , WytWy*m, WztWz*m
% 2- GtWdtWdG*m
% 3- The smallest term WstWs*m
% 4- The right-hand-side G'*WdtWd*d
% r=(A*(x0)-RHS);

%wd = std(d);

LHS = compLHS_v3(m0, G, Wr, wd, nx, ny, nz) * beta;
RHS = G' * d;

r = (LHS - RHS);
p = -r;

rold = r' * r;
%rnorm = norm(r);
%dnorm = norm(RHS);
count = 1;
%error = rnorm / dnorm;
%misfit(count) = norm(r)^2;

while (count) < 5%error>=10e-4
tic
    Ap = compLHS_v3(p, G, Wr, wd, nx, ny, nz) * beta;
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
toc

end
