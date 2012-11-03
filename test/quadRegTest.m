clear all; close all;

%%

% Program a quadratic regularization inversion to check forward operator
% accuracy against analytic solutions

addpath ../utils

N = 6:2:40;
N = 14;

%% Graded Mesh (dumb octree)
curve = 6;
weight = 0.1; % Percent finer mesh at finest scale to coarsest.

% Project cell walls to centres
centre = @(x) x(1:end - 1) + 0.5*diff(x);
% Get function for creating meshes with given arguments. See meshfunc for
% details.
gmesh = meshfunc(curve, weight);
% gmesh = meshgui();
[z, dz] = gmesh(linspace(5, 6, N - 1));
[x, dx] = gmesh(linspace(5, 6, N));
[y, dy] = gmesh(linspace(5, 6, N + 1));
    
zc = centre(z);
xc = centre(x);
yc = centre(y);

plotMesh(z, x)
   
nz = length(zc);
nx = length(xc);
ny = length(yc);
nnn = nz * nx * ny;
   
Nnn = nnn;
      
[ObsX, ObsY] = meshgrid(5 : 1/10 : 6, 5 : 1/10 : 6);
ObsZ = zeros(size(ObsX)) + 0.5;
nobs = numel(ObsX(:));   

m = zeros(nx,ny,nz);
m(:,:,6) = 1;
    
% Build the dz * dx * dy cell volume vector
h = dz' * dx;
h = h(:) * dy;
h = h(:)';
    
%%
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
G2 = Z .* H .* 1./R;
surf(reshape(G2*m(:),nx,ny))