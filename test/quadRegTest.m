clear all; close all;

%%

% Program a quadratic regularization inversion to check forward operator
% accuracy against analytic solutions

% Project cell walls to centres
centre = @(x) x(1:end - 1) + 0.5*diff(x);




x = linspace(0,1,22);
y = linspace(0,1,20);
z = linspace(2,3,18);

zc = centre(z);
xc = centre(x);
yc = centre(y);
    
    nz = length(zc);
    nx = length(xc);
    ny = length(yc);
    nnn = nz * nx * ny;   

   

dx = x(2) - x(1);
dy = y(2) - y(1);
dz = z(2) - z(1);

h = ones(1, nnn) * dx*dy*dz;

%[Z X Y] = meshgrid(z, x, y);

nobsx = 20;
[ObsX, ObsY] = meshgrid(linspace(0, 1, nobsx),...
    linspace(0,1,nobsx));

ObsZ = zeros(size(ObsX)) + 0.5;

nobs = numel(ObsX(:));      
% When I ordered the matrices to accomodate the matlab vectorize order (
% for when we multiply m(:).  the kronecker products required to build the
% operator showed a slightly different symmetry than when BenP derived the
% same problem, so the codes look a little different.


% Geometric part of the Gravity kernel
% According to the UBC Documentation, the gravity kernel involves source
% location minus observation location, whereas the magnetic kernel
% has the order reversed.  Below is the magnetic.

% R = ((((kron3(ez,exy',obsX',ez')-kron3(ez,exy,obsX,ez)).^2)...
% +(kron3(ez,exy',ez',obsY')-kron3(ez,exy,ez,obsY)).^2)...
% +(kron3(1,exy',ez',obsZ)-kron3(1,exy,z,exy)).^2).^(-1/2);
% exy = ones(1,nx*ny)
% ez = ones(1,nz)
% Here we have the gravity case.
% G = (dv*(kron3(1,exy',z,exy)-kron3(1,exy,ez',obsZ))...
% .*((((kron3(ez,exy',obsX,ez)-kron3(ez,exy,obsX',ez')).^2)...
% +(kron3(ez,exy',ez,obsY)-kron3(ez,exy,ez',obsY')).^2)...
% +(kron3(1,exy',z,exy)-kron3(1,exy,ez',obsZ)).^2).^(-3/2));



e = @ (n) ones(1, n);
% Scale up dimensions to 3D
kron3 = @(a, b, c)  kron(a, kron(b, c));
% Kron observations to size G
kronobs = @(obs) kron( e(nnn), obs(:) );
% Kron dimensions to size G and subtract obs distance
H = kron( e(nobs)', h);
Z = kron( e(nobs)', kron3( e(ny), e(nx), zc)) - kronobs(ObsZ);
X = kron( e(nobs)', kron3( e(ny), xc, e(nz))) - kronobs(ObsX);
Y = kron( e(nobs)', kron3( yc, e(nx), e(nz))) - kronobs(ObsY);
R = (X.^2 + Y.^2 + Z.^2).^(3/2);
G = Z .* H .* 1./R;

% G = dv*(kron(e(nx*ny)',Z(:)')-zeros(nx*ny,nx*ny*nz))...
% .*((((kron(e(nx*ny)',X(:)')-kron(e(nx*ny*nz),kron(x',e(ny)'))).^2)...
% +((kron(e(nx*ny)',Y(:)')-kron(e(nx*ny*nz),kron(e(nx)',y'))).^2)...
% +((kron(e(nx*ny)',Z(:)')-zeros(nx*ny,nx*ny*nz)).^2)).^(-3/2));

m = ones(nz, nx, ny);
%noise = 5 * randn(nz, nx, ny);
m(round(0.5*nz),...
    round(0.5*nx) : round(0.5*nx) + 1,...
    round(0.5*ny) - 1 : round(0.5*ny) + 1 ) = 4;

%m = m + noise;

% m(1:2,1:2,15) = 2000;
d = G*m(:);
d = d + 0.1 * mean(d) * randn(length(d), 1);
d = reshape(d, size(ObsX) );
figure()
imagesc(xc,yc,d)
title(sprintf('max data = %f', max(d(:))))
axis square
% m(1:2,1:2,15) = 2000;
D = @(n) spdiags([-e(n )', e(n )'], [-1 0], n +1 , n);

Iz = speye(nz);
Ix = speye(nx);
Iy = speye(ny);
Dz = kron(Iy, kron(Ix, D(nz)));
Dx = kron(Iy, kron(D(nx), Iz));
Dy = kron(D(ny), kron(Ix, Iz));

GRAD = [Dz; Dx; Dy];

Im = spdiags(e(nnn)' * h, 0, nnn, nnn);

beta = 1e-7;

A = G' * G + beta * (GRAD' * GRAD) + beta * Im;
minv = A\(G' * d(:));

minv = reshape(minv, nz, nx, ny);

figure()
for zi = 1:length(zc)
    subplot(round(sqrt(nz)), round(sqrt(nz)), zi)
    imagesc(squeeze(minv( zi ,:, :)))
    title(sprintf('depth %i', zi))
end