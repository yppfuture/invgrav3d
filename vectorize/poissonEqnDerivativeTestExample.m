%3D Poisson Driver
clear all
close all

%% Set constants

nn =  8;
L = 1;

n1 = nn; % Number of cells in x1 direction
n2 = nn; % Number of cells in x2 direction
n3 = nn; % Number of cells in x3 direction

nnn = n1 * n2 * n3;

h1 = L/n1; % cell length in x1 direction
h2 = L/n2;
h3 = L/n3;

Lx1 = n1*h1;
Lx2 = n2*h2;
Lx3 = n3*h3;

%% Set up cells, mesh, Operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MESH
[x,y,z] = ndgrid(0:h1:Lx1,0:h2:Lx2,0:h3:Lx3); % Cell nodes
[xc,yc,zc] = ndgrid(h1/2:h1:Lx1-h1/2, h2/2:h2:Lx2-h2/2, h3/2:h3:Lx3-h3/2); % Cell centres
[xdx,ydx,zdx] = ndgrid(0:h1:Lx1, h2/2:h2:Lx2-h2/2, h3/2:h3:Lx3-h3/2); %Staggered in x1 cell wall x2 x3
[xdy,ydy,zdy] = ndgrid(h1/2:h1:Lx1-h1/2, 0:h2:Lx2, h3/2:h3:Lx3-h3/2); %Staggered in x2 cell wall x1 x3
[xdz,ydz,zdz] = ndgrid(h1/2:h1:Lx1-h1/2, h2/2:h2:Lx2-h2/2, 0:h3:Lx3); %Staggered in x3 cell wall x1 x2

%%%%%%%%%%%%%%%%%% OPERATORS
% derivative function (calculating from nodes to centre) <DIV>
ddxc = @(m,k) 1/k*spdiags([-ones(m+1,1) ones(m+1,1)],[0,1],m,m+1); 
% Average function (Go from centres to walls (will also need ghost))
av = @(m) 0.5*spdiags([ones(m,1) ones(m,1)],[-1,0],m+1,m); 
%%%%%%% 

I1   = speye(n1);     % Create Identities of Appropriate size
I2   = speye(n2);     % Create Identities of Appropriate size
I3   = speye(n3);     % Create Identities of Appropriate size

Dc1  = ddxc(n1, h1);  % Create 1D Operators 
Dc2  = ddxc(n2, h2);  % Create 1D Operators 
Dc3  = ddxc(n3, h3);  % Create 1D Operators 

Av1  = av(n1);
Av2  = av(n2);
Av3  = av(n3);

% Boundary Conditions
Av1([1, end]) = 2 * [1,1];
Av2([1, end]) = 2 * [1,1];
Av3([1, end]) = 2 * [1,1];

%%% Ramp up to 2D then to 3D
%% 3D
DC1 = kron(I3, kron(I2, Dc1));
DC2 = kron(I3, kron(Dc2, I1));
DC3 = kron(Dc3, kron(I2, I1));

A1 = kron(I3, kron(I2, Av1));
A2 = kron(I3, kron(Av2, I1));
A3 = kron(Av3, kron(I2, I1));
AV = [A1; A2; A3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Derivative Test
% Scale up operators to work in block diag dimension
D  = [DC1, DC2, DC3];
% Create a bunch of random vectors for the test.
u = randn(size(D',2),1);
m = randn(size(u));
q = randn(size(m));
v = randn(size(m));

% Sparse diag of 1/m^2
diagm = @(m) spdiags(1./m.^2, 0, nnn, nnn) ;

% Inverse of Avg*1/m diagonalized ie diag( 1/(Avg*(1/m)) )
Ainv = @(A, m) ( spdiags(1./(A*(1./m)), 0, size(A, 1), size(A, 1)) );
% Big block diag of the inverse averaged m operator above
diagAvgm = @(m) (blkdiag(Ainv(A1,m), Ainv(A2, m), Ainv(A3, m)));

% Forward Simulation Operator C
C = @(u,m) (D * diagAvgm(m) * D' * u - q);
% Derivate Functions
dCdu = @(u,m)( D * diagAvgm(m) * D' );
dCdm = @(u,m)( D * diag(D'*u) * (diagAvgm(m).^2) * AV * diagm(m));

% Wrappers to turn funcs into 1 arg funcs to pass into derivTest func
Cu = @(u) ( C(u, m) );
Cm = @(m) ( C(u, m) );

% Perform deriv test for dCdu
fprintf('_____________________________________\n')
fprintf('Calculating Derivative test for dC/du\n')
dCu = @(u) ( dCdu(u, m) );
derivTest(Cu, dCu, length(u));
% Perform deriv test for dCdm
fprintf('_____________________________________\n')
fprintf('Calculating Derivative test for dC/dm\n')
dCm = @(m) ( dCdm(u, m) );
derivTest(Cm, dCm, length(u));

