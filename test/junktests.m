% Graded mesh
clear all
close all

addpath ../vectorize
addpath ../utils
kron3 = @(a, b, c)  kron(a, kron(b, c));
e = @ (n) ones(n, 1);
%% gmesh tests
%{
n = 50;

curve = 4;
w = 0.1;
% Gmesh builds a fat parabola with minimum dictated by wieght
% The w/(1-w) is necessary to that the scaling is preserved.
gmesh = @(x) ((x - (x(end) - x(1))/2 - x(1)).^curve + (w / (1-w) ) * ((x(end) - x(1))/2).^curve);
% getdx turns gmesh into an x and dx
getdx = @(x) (x(end) - x(1)) * gmesh(x) / sum(gmesh(x));
%getdx = @(x) gmesh(x);
getx = @(x) cumsum([x(1), getdx(x)]);
% Project cell walls to centres
centre = @(x) x(1:end - 1) + 0.5*diff(x);
% Get dims performs a few function calls to make things easier.
getdims = @(x) deal( getx(centre(x)), getdx(centre(x)) );



gd = meshfunc(4, 0.1, true);

x0 = linspace(2, 6 , n);

[x, dx] = gd( x0 );
xc = centre(x);
[x2, dx2] = getdims( x0 );

norm(x - x2)




fprintf('Ratio of dx(centre) / dx(end) = %f\n', dx(round(n/2))/dx(end))

%plot(xc, [dx], 'r*')
%line([x; x], [ones(n,1)' * min(dx) ; ones(n,1)' *max(dx)],'Color',[.8 .8 .8])
%xlim( [min(x) - 0.01*min(x), max(x) + 0.01*min(x)] )
%ylim( [min(dx), max(dx) ] )
%title('Graded mesh')
%legend('Calculated dx','mesh line \rho')

%}

%% Isosurface model tests
% More Realistic Models

nx = 15;
ny = 15;
nz = 15;

% Try doing this on the graded mesh...

rng('shuffle', 'twister')
data = rand(nz-10, nx-4, ny-4);
data = smooth3(data,'box', 1);
m = 0.5*max(max(max(data)));
data(data < m) = 0;
data = 3*data;
rng(5,'twister')
M = 0.5*rand(nz, nx, ny);
M(4:8, 3:end-2, 3:end-2) = data;
M = smooth3(M,'box', 3);

m = 0.7*max(max(max(M)));
% Too keep whole grid in scope... prob a better way.
M(1,1,1) = 2;
M(end,end,end) = 2;

figure()
patch(isocaps(M, m),'FaceColor','interp','EdgeColor','none');
p1 = patch(isosurface(M, m),'FaceColor','blue','EdgeColor','none');
isonormals(M,p1)
%zlim([0,nz])
%xlim([0,nx])
%ylim([0,ny])
view(3); 
axis vis3d tight
camlight left; 
lighting phong
axis square
grid on

%figure()

%imagesc(squeeze(dw(:,:,1)))
%xlabel('depth')
%ylabel('x dim')