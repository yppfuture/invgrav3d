clear all
close all
% Set up dimensions

N = 6:2:10;
%N = 10;
% Outer Test loop
for tt = 1:length(N)
    x = linspace(5, 6, N(tt))';
    y = linspace(5, 6, N(tt) + 1)';
    z = linspace(5, 6, N(tt) - 1)';
    
    n1 = length(z);
    n2 = length(x);
    n3 = length(y);
    nnn = n1 * n2 * n3;
    
    Nnn(tt) = nnn;
    
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    dz = z(2) - z(1);
    
    [ObsX, ObsY] = meshgrid(5 : 1/10 : 6, 5 : 1/10 : 6);
    ObsZ = zeros(size(ObsX)) + 0.5;
    nobs = numel(ObsX(:));
    
    h = dx*dy*dz;
    
    
    %% Loop Code
    tic
    g = zeros(n1, n2, n3);
    G = zeros(nnn, nobs);
    count = 1;
    for kk = 1:nobs
        for iy = 1 : n3
            for ix = 1 : n2
                for iz = 1 : n1
                    g(iz, ix, iy) =((z(iz) - ObsZ(kk) )*h)*1/(((x(ix)-ObsX(kk))^2 + (y(iy) - ObsY(kk))^2 + (z(iz) - ObsZ(kk))^2).^(3/2));
                end
            end
        end
        G(:,count) = g(:);%reshape(g,[nnn,1]);
        count = count +1;
    end
    loopt(tt) = toc; %#ok<*SAGROW>
    
    %% Vectorized Code
    tic
    e = @ (n) ones(n, 1);
    % Scale up dimensions to 3D
    kron3 = @(a, b, c)  kron(a, kron(b, c));
    % Kron observations to size G
    kronobs = @(obs) kron( e(nnn), obs(:)' );
    % Kron dimensions to size G
    Z = (kron( e(nobs)', kron3( e(n3), e(n2), z)) - kronobs(ObsZ));
    X = (kron( e(nobs)', kron3( e(n3), x, e(n1))) - kronobs(ObsX));
    Y = (kron( e(nobs)', kron3( y, e(n2), e(n1))) - kronobs(ObsY));
    R = (X.^2 + Y.^2 + Z.^2).^(3/2);
    G2 = Z * h .* 1./R;
    
    vect(tt) = toc;
     
end

figure()
fprintf('Norm of (Loop code - Vec code) is %f\n', norm(G - G2))
plot(Nnn, loopt, Nnn, vect)
legend('Loop time', 'Vectorized time')
title('Computation Time Comparison of \nVectorized vs Loop based Forward Operator Code')
ylabel('Time [seconds]')
xlabel('Number of cells')