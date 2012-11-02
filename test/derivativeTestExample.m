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

%% Operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% derivative function (calculating from nodes to centre) <DIV>
ddxc = @(m,k) 1/k*spdiags([-ones(m+1,1) ones(m+1,1)],[0,1],m,m+1); 
% Average function (Go from centres to walls (will also need ghost))
av = @(m) 0.5*spdiags([ones(m,1) ones(m,1)],[-1,0],m+1,m); 

Dc1  =  ddxc(n1, h1);  % Create 1D Operators 
Dc2  = ddxc(n2, h2);  % Create 1D Operators 
Dc3  = ddxc(n3, h3);  % Create 1D Operators 

Av1  = av(n1);
Av2  = av(n2);
Av3  = av(n3);
 
% Boundary Conditions
Av1([1, end]) = 2 * [1,1];
Av2([1, end]) = 2 * [1,1];
Av3([1, end]) = 2 * [1,1];

%% Scale to 3D
% Convenience funcs
I = @(n) speye(n);
kron3 = @(a, b, c)  kron(a, kron(b, c));

% 3d grad/divergence operators
DC1 = kron3(I(n3), I(n2), Dc1);
DC2 = kron3(I(n3), Dc2, I(n1));
DC3 = kron3(Dc3, I(n2), I(n1));

% 3d averaging operators
A1 = kron3(I(n3), I(n2), Av1);
A2 = kron3(I(n3), Av2, I(n1));
A3 = kron3(Av3, I(n2), I(n1));

AV = [A1; A2; A3];
D  = [DC1, DC2, DC3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Derivative Test
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

