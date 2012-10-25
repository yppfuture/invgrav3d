%%
clear all; close all
%%

%%
% The main script to perform simplified vectorized gravity inversion
%%

%% Parameter definition

% Let's start in 1D

dx = 0.1          % x - vector defined by length Lx 
Lx = 100;         % and increment length dx
nx = Lx/dx;
x = 0:dx:Lx;

% 3D

y = x; dy = dx;   % same vectors and same increments
z = x; dz = dx;

cube = ndgrid(x,y,z)
%%