clear all; close all;

%%

% Program a quadratic regularization inversion to check forward operator
% accuracy against analytic solutions

Lx = 10;
Ly = 15;
Lz = 20;


nx = 20;
ny = 10;
nz = 15;

x = linspace(1,Lx,nx);
y = linspace(1,Ly,ny);
z = linspace(1,Lz,nz);

dx = Lx/nx;
dy = Ly/ny;
dz = Lz/nz;

dv = dx*dy*dz;

[X Y Z] = meshgrid(x,y,z);

% When I ordered the matrices to accomodate the matlab vectorize order (
% for when we multiply m(:).  the kronecker products required to build the
% operator showed a slightly different symmetry than when BenP derived the
% same problem, so the codes look a little different.

kron3 = @(a,b,c,d) kron(a,kron(b,kron(c,d)));
e = @(n) ones(1,n);

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

G = dv*(kron(e(nx*ny)',Z(:)')-zeros(nx*ny,nx*ny*nz))...
.*((((kron(e(nx*ny)',X(:)')-kron(e(nx*ny*nz),kron(x',e(ny)'))).^2)...
+((kron(e(nx*ny)',Y(:)')-kron(e(nx*ny*nz),kron(e(nx)',y'))).^2)...
+((kron(e(nx*ny)',Z(:)')-zeros(nx*ny,nx*ny*nz)).^2)).^(-3/2));

m = zeros(nx,ny,nz);
m(5:5,5) = 2000;
% m(1:2,1:2,15) = 2000;

d = reshape(G*m(:),nx,ny);
figure()
imagesc(x,y,d)


% e = ones(nx+1,1)
% D = spdiags([-e e],[0 1],nx,nx)

