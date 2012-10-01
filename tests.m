% Symmetric matrix test
clear all
close all


n = 5;
m = 10;
A = randi(10,[n,m]);
x = randi(10, [m,1]);
full(A')

B = A' * A;

full(B)