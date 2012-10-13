%Inverting for compact model - Iterative method
%Dominique Fournier
%October 15, 2011

clear all
close all
addpath data
addpath functions
addpath ../matlab/imStacks

load data;
load kernel;
load model;
load Wr;

%Extract cell dimensions from UBC mesh file
[dx, dy, dz] = getMesh( importdata('UBC_mesh.msh', ' ', 0) );

nx = length(dx);
ny = length(dy);
nz = length(dz);

Lx = min(dx);
Ly = min(dy);
Lz = min(dz);

% background=2.7; %Background density

mcell = nx * ny * nz;
mdata = length(data);
% target=0.5;     %Expected density contrast

%Modify data and model to substract background density
% model=zeros(size(rho,1)*size(rho,2),1);
% model(:)=rho(:)-background;

m0 = zeros(mcell, 1);
% m0(:)=background;

%% Compute compact model using iterative method
% alphaS=1;
% alphax=alphaS*Lx^2;
% alphaz=alphaS*Lz^2;

alphaS = 1 / Lx.^2;
alphax = 1.0;
alphay = 1.0;
alphaz = 1.0;

% Define coefficients

beta=1e-1; % Trade-off parameter
wd=std(data).^(1/2);
invm = conjGrad(m0, G, data, Wr, wd, nx, ny, nz, beta);
save('data/model_out.dat','-ascii','invm')

misfit=sqrt(sum((G*invm-mdata).^2));
% Slicer(reshape(invm, nz, nx, ny))

%Slicer(reshape(m, nz, nx, ny))