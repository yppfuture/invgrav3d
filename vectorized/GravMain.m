%%

clear all; close all

%%

% The main script to perform simplified vectorized gravity inversion

%% Parameter definition

% Let's start in 1D

dx = 1;          % x - vector defined by length Lx 
Lx = 10;         % and increment length dx
nx = Lx/dx;
x = 1:dx:Lx;

% 3D

y = x; dy = dx; Ly = Lx; ny = nx;   % same vectors and same increments
z = x; dz = dx; Lz = Lx; nz = nx;

dV = dx*dy*dz; % Volume element


n1 = nz; n2 = nx; n3 = ny;

nnn = n1 * n2 * n3;

[ObsX ObsY] = meshgrid(x,y);
ObsZ = zeros(size(ObsX));
nobs = numel(ObsX(:));

% Func to creat n-dimensional ones vector
e = @(n) ones(1 ,n);

% Func to scale up dimensions to 3D
kron3 = @(a, b, c)  kron(a, kron(b, c));

% Func to scale observations to size G
kronobs = @(obs) kron( e(nnn), obs(:) );

% Produce a G matrix 
Z = (kron( e(nobs)', kron3( e(n3), e(n2), z)) - kronobs(ObsZ));
X = (kron( e(nobs)', kron3( e(n3), x, e(n1))) - kronobs(ObsX));
Y = (kron( e(nobs)', kron3( y, e(n2), e(n1))) - kronobs(ObsY));
R = (X.^2 + Y.^2 + Z.^2).^(3/2);
G = (6.67e-11)*(Z * dV .* 1./R);

%% Test the operator for the nx = ny = nz = 10 case
 
% rho = zeros(nx,ny,nz);
% rho(:,5:6,5) = 1;      % add some buried feature of anomalous density
% 
% d = G*rho(:);  % G on rho (model) = data
% figure(1)
% contourf(reshape(d,nx,ny))
% 
% %% Corrupt the signal
% 
% dNoisy = (0.05*randn(size(d))+(0.01*min(d))) +d;
% figure(2)
% contourf(reshape(dNoisy,nx,ny))

%%  Getting there

% Build the gradient operator *** - rework this for z=3 ***

d1 = spdiags([-e(nx)' e(nx)'],[0 1],nx,nx);
d2 = spdiags([-e(ny)' e(ny)'],[0 1],ny,ny);
d3 = spdiags([-e(nz)' e(nz)'],[0 1],nz,nz);

D1 = kron3(d1,speye(ny),speye(nz)); % x dimension
D2 = kron3(speye(nx),d2,speye(nz)); % y dimension
D3 = kron3(speye(nx),speye(ny),d3); % z dimension

DDD = [D1;D2;D3]; % 3Dify

av1 = spdiags([e(nx)' e(nx)'],[0 1],nx,nx);
av2 = spdiags([e(ny)' e(ny)'],[0 1],ny,ny);
av3 = spdiags([e(nz)' e(nz)'],[0 1],nz,nz);

AV1 = kron3(av1,speye(ny),speye(nz));
AV2 = kron3(speye(nx),av2,speye(nz));
AV3 = kron3(speye(nx),speye(ny),av3);

AV = [AV1;AV2;AV3];

% Build an A for Am=b to be used in conjugate gradient minimization

beta = 1; % set this to one for now
A = G'*G + beta*DDD'*diag(AV*e(nnn)')*DDD;

% Build b for Am=b

% b = G'd









