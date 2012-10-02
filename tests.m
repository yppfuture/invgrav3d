% Symmetric matrix test
clear all
close all

addpath data
load kernel
load model

N = size(G, 2);
M = size(G, 1);

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

