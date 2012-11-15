function GRAD = gradientOp(n1, n2, n3, dx1, dx2, dx3)
% GRADIENTOP produces a 3 dimensional gradient operator
% Note that GRAD' is divergence.
% Grad takes cell nodes and computes derivatives on faces (vectors)
% GRAD' or DIV takes values on walls and computes nodes (scalars)

e = @ (n) ones(n, 1);
kron3 = @(a, b, c)  kron(a, kron(b, c));
D = @(n,k) spdiags([-e(n) .* k, e(n) .* 1/k], [-1 0], n + 1 , n);

I1 = speye(n1);
I2 = speye(n2);
I3 = speye(n3);

GRAD = [
    kron(I3, kron(I2, D(n1, dx1)))
    kron(I3, kron(D(n2, dx2), I1))
    kron(D(n3, dx3), kron(I2, I1))
    ];

end
