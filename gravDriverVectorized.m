clear all; close all;


addpath test
addpath vectFuncs
addpath functions
%%

% Program a quadratic regularization inversion to check forward operator
% accuracy against analytic solutions

% Project cell walls to centres
centre = @(x) x(1:end - 1) + 0.5*diff(x);

x = linspace(0,1,10);
y = linspace(0,1,10);
z = linspace(2,3,10);

x = x(:); y = y(:); z = z(:);

zc = centre(z);
xc = centre(x);
yc = centre(y);

nz = length(zc);
nx = length(xc);
ny = length(yc);

nnn = nz * nx * ny;

dx = diff(x);
dy = diff(y);
dz = diff(z);

dxc = diff(xc);
dyc = diff(yc);
dzc = diff(zc);

%[Z X Y] = meshgrid(z, x, y);

nobsx = 20;
[ObsX, ObsY] = meshgrid(linspace(0, 1, nobsx),...
    linspace(0,1,nobsx));

ObsZ = zeros(size(ObsX)) + 0.5;

nobs = numel(ObsX(:));

e = @ (n) ones(n, 1);
% Scale up dimensions to 3D
kron3 = @(a, b, c)  kron(a, kron(b, c));
% Kron observations to size G
kronobs = @(obs) kron( e(nnn)', obs(:) );

% Kron up dz dy dz and get cell volume h
h = kron3(e(ny), e(nx), dz(:)) .* ...
    kron3(e(ny), dx(:), e(nz)) .* ...
    kron3(dy(:), e(nx), e(nz)) ;

% Kron dimensions to size G and subtract obs distance
H = kron( e(nobs), h');
Z = kron( e(nobs), kron3( e(ny)', e(nx)', zc')) - kronobs(ObsZ);
X = kron( e(nobs), kron3( e(ny)', xc', e(nz)')) - kronobs(ObsX);
Y = kron( e(nobs), kron3( yc', e(nx)', e(nz)')) - kronobs(ObsY);
R = (X.^2 + Y.^2 + Z.^2).^(3/2);
G = Z .* H .* 1./R;



m = ones(nz, nx, ny);
m(round(0.5*nz),...
    round(0.5*nx) : round(0.5*nx) + 1,...
    round(0.5*ny) : round(0.5*ny) + 1 ) = 4;



d = G*m(:);
d = d + 0.05 * mean(d) * randn(length(d), 1);
dcube = reshape(d, size(ObsX) );
figure()
imagesc(xc,yc,dcube)
title(sprintf('max data = %f', max(d(:))))
axis square


GRAD = gradientOp(nz, nx, ny, dzc, dxc, dyc);

Im = spdiags(e(nnn) * h(1), 0, nnn, nnn);

beta = 1e-10;

A = G' * G + beta * (GRAD' * GRAD) + beta * Im;
minv = A\(G' * d(:));

minv = reshape(minv, nz, nx, ny);

figure()
for zi = 1:length(zc)
    subplot(round(sqrt(nz)), round(sqrt(nz)), zi)
    imagesc(squeeze(minv( zi ,:, :)))
    title(sprintf('depth %i', zi))
end



%% Derivative Testing
% Phi = 1/2 ||Gm-d||^2 + 1/2(|grad*m|^2 + m'm*h)
% dPhi/dm = G'Gm - G'd + grad'grad*m + m*h
m = m(:);
d = d(:);

mr = randn(size(G,2),1);
Phi = @(m) 1/2 * (G*m - d)' * (G*m -d) + ...
    1/2 * ( (GRAD * m)' * (GRAD * m)) ;
    

dPhidm = @(m) ( (G' * G) * m) - (G' * d) + ...
    (GRAD' * GRAD) * m;


derivTest(Phi, dPhidm, length(mr))




