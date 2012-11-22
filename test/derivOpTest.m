clear all; close all;

addpath ../vectFuncs

centre = @(x) x(1:end - 1) + 0.5*diff(x);

x1 = linspace(0,10,100);
xc = centre(x1);

dx1 = diff(xc); 
dx2 = dx1; 
dx3 = dx1;

[XC1 XC2 XC3] = meshgrid(xc);
[X1 X2 X3]    = meshgrid(x1);


n1 = length(xc); 
n2 = n1; 
n3 = n1;
fc = sin(pi*XC1).*sin(pi*XC2).*sin(pi*XC3);
f  = sin(pi*X1).*sin(pi*X2).*sin(pi*X3);
% GRAD = gradientOp(n1, n2, n3, dx1(:), dx2(:), dx3(:));
dx1 = dx1(:); dx2 = dx2(:); dx3 = dx3(:);
dx1 = [dx1(1); dx1; dx1(end)];
dx2 = [dx2(1); dx2; dx2(end)];
dx3 = [dx3(1); dx3; dx3(end)];

e = @(n) ones(n, 1);
kron3 = @(a, b, c)  kron(a, kron(b, c));
D1 = kron( dx1, e(n1)') .* spdiags([-e(n1) , e(n1) ], [-1 0], n1 + 1 , n1);
D2 = kron( dx2, e(n2)') .* spdiags([-e(n2) , e(n2) ], [-1 0], n2 + 1 , n2);
D3 = kron( dx3, e(n3)') .* spdiags([-e(n3) , e(n3) ], [-1 0], n3 + 1 , n3);

I1 = speye(n1);
I2 = speye(n2);
I3 = speye(n3);

GRAD = [
    kron(I3, kron(I2, D1))
    kron(I3, kron(D2, I1))
    kron(D3, kron(I2, I1))
    ];


% 
% GRADf = [pi*cos(pi*X1).*sin(pi*X1).*sin(pi*X1);
%          pi*sin(pi*X1).*cos(pi*X2).*sin(pi*X3);
%          pi*sin(pi*X1).*sin(pi*X2).*cos(pi*X3)];

LAPf = -3*pi*pi*sin(pi*XC1).*sin(pi*XC2).*sin(pi*XC3);

LAPresid = norm( LAPf(:) - (GRAD'*GRAD)*fc(:) )
% GRADresid = norm( GRADf(:)-GRAD*fc(:) );


