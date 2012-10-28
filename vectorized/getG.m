function [G] = getG(n2,n3,n1,x,y,z,nobs,ObsX,ObsY,ObsZ)

nnn = n2*n3
e = @ (n) ones(n, 1);
% Scale up dimensions to 3D
kron3 = @(a, b, c)  kron(a, kron(b, c));
% Kron observations to size G
kronobs = @(obs) kron( e(nnn), obs(:)' );
% Kron dimensions to size G
Z = (kron( e(nobs)', kron3( e(n3), e(n2), z)) - kronobs(ObsZ));
X = (kron( e(nobs)', kron3( e(n3), x, e(n1))) - kronobs(ObsX));
Y = (kron( e(nobs)', kron3( y, e(n2), e(n1))) - kronobs(ObsY));
R = (X.^2 + Y.^2 + Z.^2).^(3/2);
G = Z * h .* 1./R;