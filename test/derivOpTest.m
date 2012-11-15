clear all; close all;

centre = @(x) x(1:end - 1) + 0.5*diff(x);
x1 = linspace(0,10,11);
xc = centre(x1);

dx1 = diff(x1); dx2 = dx1; dx3 = dx1;
[X1 X2 X3] = meshgrid(xc);
n1 = length(x1); n2 = n1; n3 = n1;
f = sin(pi*X1).*sin(pi*X2).*sin(pi*X3);

GRAD = gradientOp(n1, n2, n3, dx1, dx2, dx3)

GRAD*f(:)


GRADf = [pi*cos(pi*X1).*sin(pi*X2).*sin(pi*X3);
         pi*sin(pi*X1).*cos(pi*X2).*sin(pi*X3);
         pi*sin(pi*X1).*sin(pi*X2).*cos(pi*X3)]

DIV = GRAD'*GRAD;


DIVf = -3*pi*pi*sin(pi*X1).*sin(pi*X2).*sin(pi*X3);

DIVresid = sum((DIVf-DIV*f(:)).^2);
GRADresid = sum((GRADf-GRAD*f(:)).^2);