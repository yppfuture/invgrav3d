% Graded mesh
clear all
close all
addpath vectorize


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

