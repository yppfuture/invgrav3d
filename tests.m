% Graded mesh
clear all
close all

n = 100;
x = linspace(2, 6 , n);

weight = 5;
% Gmesh builds a fat parabola with minimum dictated by wieght
gmesh = @(x) ((x - (x(end) - x(1))/2 - x(1)).^4 + (x(end) - x(1)) / (length(x) / weight));
% getdx turns gmesh into an x and dx
getdx = @(x) (x(end) - x(1)) * gmesh(x) / sum(gmesh(x));
getx = @(x) cumsum([x(1), getdx(x)]);

dx = getdx(x);
x = getx(x);

plot(dx)