% Symmetric matrix test
clear all
close all

addpath data
addpath functions
%load kernel
%
%N = size(G, 2);
%M = size(G, 1);
% 
% m = randn(N,1);
% 


%{
tic
GtG_rowii = zeros(1, N);
GtGm1 = zeros(N, 1);
for ii = 1 : N
    for jj = 1 : N
        GtG_rowii(jj) = 0;
        for kk = 1 : M
            GtG_rowii(jj) = GtG_rowii(jj) + G(kk, ii) * G(kk, jj);
        end
        GtGm1(ii) = GtGm1(ii) + GtG_rowii(jj) * m(jj);
    end
    
end
toc

tic
GtG_rowii = zeros(1, N);
GtGm2 = zeros(N, 1);
for ii = 1 : N % Main Loop
    for jj = ii : N  % Note loop runs from ii -> N
        GtG_rowii(jj) = 0; % Clear previous results
        for kk = 1 : M % Dot product of G columns
            GtG_rowii(jj) = GtG_rowii(jj) + G(kk, ii) * G(kk, jj);
        end
        GtGm2(ii) = GtGm2(ii) + GtG_rowii(jj) * m(jj); % Compute top half of triangle
    end
    if ii + 1 <= N
        for ll = ii + 1 : N 
            GtGm2(ll) = GtGm2(ll) + GtG_rowii(ll) * m(ii); % Compute lower half of triangle
        end
    end
end
toc

err1 = sqrt( (GtGm2 - GtGm1)' * (GtGm2 - GtGm1) ) ;
%}
X0 = 0;
Y0 = 0;
Z0 = 0;
Nx = 10:19;
Ny = Nx;
Nz = 4:13;
dz = 1;
dx = 1;
dy = 1;
% Speed profile for the fun of it:
for zz = 1 : 10

nx = Nx(zz);
ny = nx;
nz = Nz(zz);
    

count = 1;
[ObsX, ObsY] = meshgrid(5:12,5:12);
ObsZ = 1;
n = size(ObsX, 2) * size(ObsY, 1);
G = zeros(n,nx*ny*nz);
for ii = 1 : n;
        G(count,:) = forwardGrav(nx, ny, nz, X0, Y0, Z0,...
            dx, dy, dz, ObsX(ii), ObsY(ii), ObsZ);
count = count + 1;
end

N = size(G, 2);
M = size(G, 1);

m = randn(N,1);

tic
[~] = compLHS2(m, G, nx, ny, nz);
cl2(zz) = toc;
tic
[~] = compLHS(m, G, nx, ny, nz);
cl(zz) = toc;
tic
[~] = getLHS(m, G, nx, ny, nz);
gl(zz) = toc;
tic
[~] = getLHS2(m, G, nx, ny, nz);
gl2(zz) = toc;

end
%%
load cdata.mat
psize = Nx .* Ny .* Nz;
bigp = 1 : 100 : 10000;
% func for plotting
g = @(x,t)( x(1) * exp( x(2) * t) );

% Linearize data -> build single diagonal linear system
% Note t = x1 * exp( x2 * psize) -> ln t = ln x1 + x2 * psize
A = kron(eye(4),[ones(length(psize),1), psize(:)]);
t = log([cl,cl2,gl,gl2]');
% Solve equations (At*A \ At*b) and delinearize x1. 
x = (A' * A) \ (A' * t);
x([1, 3, 5, 7]) = exp(x([1, 3, 5, 7]));

figure(1)
subplot(3,1,1)
loglog(psize,cl2, 'b--', psize, cl, 'b', psize, gl, 'r--', psize, gl2, 'r')
title('LogLog algorithm time versus increasing problem size')
legend('Matlab w/o optz','Matlab w/ optz','C w/o optz', 'C w/ optz',...
    'Location','Best')
xlabel('Problem Size [nx*ny*nz]')
ylabel('Time [s]')
xlim([psize(1),psize(end)])

subplot(3,1,2)
plot(psize,cl2, 'b*', psize, cl, 'b^', psize, gl, 'r*', psize, gl2, 'r^')
title('Algorithm time versus increasing problem size')
legend('data: Matlab w/o optz','data: Matlab w/ optz',...
    'data: C w/o optz', 'data: C w/ optz',...
    'Location','Best')
xlabel('Problem Size [nx*ny*nz]')
ylabel('Time [s]')

subplot(3,1,3)
plot(bigp, g([x(1),x(2)], bigp), bigp, g([x(3),x(4)], bigp), ...
    bigp, g([x(5),x(6)], bigp), bigp, g([x(7),x(8)], bigp))
title('Projections for computational time for larger problem size')
legend('proj: Matlab w/o optz','proj: Matlab w/ optz',...
    'proj: C w/o optz', 'proj: C w/ optz',...
    'Location','Best')


%err1 = sqrt( (LHS - LHS2)' * (LHS - LHS2) ) ;

