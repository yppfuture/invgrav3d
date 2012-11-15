clear all; close all;

%%

% Program a quadratic regularization inversion to check forward operator
% accuracy against analytic solutions

%%  Set the problem size, dimensions

% Project cell walls to centres
centre = @(x) x(1:end - 1) + 0.5*diff(x);

% source position vectors on the nodes
x = linspace(0,1,22);
y = linspace(0,1,20);
z = linspace(2,3,18);

% source position vectors in centers
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
h = ones(1, nnn) * dx*dy*dz; % form volume element vector

% observation matrices
nobsx = 20;
[ObsX, ObsY] = meshgrid(linspace(0, 1, nobsx),linspace(0,1,nobsx));
ObsZ = zeros(size(ObsX)) + 0.5; % ensure observations are above model
nobs = numel(ObsX(:));          

%% Form the Forward kernel matrix

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

%% Build a model, Test an inversion
m = ones(nz, nx, ny);
%noise = 5 * randn(nz, nx, ny);
m(round(0.5*nz),...
    round(0.5*nx) : round(0.5*nx) + 1,...
    round(0.5*ny) - 1 : round(0.5*ny) + 1 ) = 4;

figure()
for zi = 1:length(zc)-1
    subplot(floor(sqrt(nz)), floor(sqrt(nz)), zi)
    imagesc(squeeze(m( zi ,:, :)))
    title(sprintf('depth %i', zi))
end

%m = m + noise;

% m(1:2,1:2,15) = 2000;
d = G*m(:);
figure()
imagesc(xc,yc,d); title('data before noise');
d = d + 0.05 * mean(d) * randn(length(d), 1) + 0.005*mean(d);
d = reshape(d, size(ObsX) );
figure()
imagesc(xc,yc,d); title('data after noise');
title(sprintf('max data = %f', max(d(:))))
axis square

% Form the Gradient operator in 3D
D = @(n) spdiags([-e(n )', e(n )'], [-1 0], n +1 , n);
Iz = speye(nz);
Ix = speye(nx);
Iy = speye(ny);                 % Boundary Conditions %
Dz = kron(Iy, kron(Ix, D(nz))); Dz(1) = 1; Dz(end) = 2;
Dx = kron(Iy, kron(D(nx), Iz)); Dx(1) = 2; Dx(end) = 2;
Dy = kron(D(ny), kron(Ix, Iz)); Dy(1) = 2; Dy(end) = 2;

GRAD = [Dz; Dx; Dy];

Im = spdiags(e(nnn)' * h, 0, nnn, nnn);

beta = 1e-7;

A = G' * G + beta * (GRAD' * GRAD) + beta * Im;
minv = A\(G' * d(:));

minv = reshape(minv, nz, nx, ny);

figure()
for zi = 1:length(zc)-1
    subplot(round(sqrt(nz)), round(sqrt(nz)), zi)
    imagesc(squeeze(minv( zi ,:, :)))
    title(sprintf('depth %i', zi))
end