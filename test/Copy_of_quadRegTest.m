clear all; close all;

%%

% Program a quadratic regularization inversion to check forward operator
% accuracy against analytic solutions

Lx = 10;
Ly = 10;
Lz = 10;

Nx = 10;
Ny = 10;
Nz = 10;

x = linspace(1,Lx,Nx);
y = linspace(1,Ly,Ny);
z = linspace(1,Lz,Nz);

dx = Lx/Nx;
dy = Ly/Ny;
dz = Lz/Nz;

dv = dx*dy*dz;

obsX = x;
obsY = y;
obsZ = zeros(Lx,Ly);
% When I ordered the matrices to accomodate the matlab vectorize order (
% for when we multiply m(:).  the kronecker products required to build the
% operator showed a slightly different symmetry than when BenP derived the
% same problem, so the codes look a little different.

kron3 = @(a,b,c,d) kron(a,kron(b',kron(c,d)));
ez = ones(1,Nz);
exy = ones(1,Nx*Ny);

% Geometric part of the Gravity kernel
% According to the UBC Documentation, the gravity kernel involves source
% location minus observation location, whereas the magnetic kernel
% has the order reversed.  Below is the magnetic.

% R = ((((kron3(ez,exy',obsX',ez')-kron3(ez,exy,obsX,ez)).^2)...
% +(kron3(ez,exy',ez',obsY')-kron3(ez,exy,ez,obsY)).^2)...
% +(kron3(1,exy',ez',obsZ)-kron3(1,exy,z,exy)).^2).^(-1/2);

% Here we have the gravity case.
G = (6.67e-11*dv*(kron3(1,exy,z,exy)-kron3(1,exy',ez',obsZ))...
.*((((kron3(ez,exy,obsX,ez)-kron3(ez,exy',obsX',ez')).^2)...
+(kron3(ez,exy,ez,obsY)-kron3(ez,exy',ez',obsY')).^2)...
+(kron3(1,exy,z,exy)-kron3(1,exy',ez',obsZ)).^2).^(-3/2));

m = zeros(Nx,Ny,Nz);
m(:,:,5) = 1;

d = reshape(G*m(:),Nx,Ny);
surf(d)
% e = ones(Nx+1,1)
% D = spdiags([-e e],[0 1],Nx,Nx)

