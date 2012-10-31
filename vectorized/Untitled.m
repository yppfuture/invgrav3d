close all; clear all
x = [1 2];
%X = meshgrid(x)
I = ones(1,2);
I4 = ones(1,4);
I8 = ones(1,8);
%kron(I4,kron(X,I4))
%L = kron(I,kron(X,I))

RHS = kron(I,kron(I4',kron(I',x)));

A = kron(I,kron(I,x))
B = kron(I4',A)

C = kron(I4,[1 1;1 1;2 2;2 2]')

B-C

