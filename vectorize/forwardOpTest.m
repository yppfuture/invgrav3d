clear all
close all
% Set up dimensions

N = 6:2:40;
%N = 10;
% Outer Test loop, for testing computational time.
for tt = 1:length(N)
    
    x = linspace(5, 6, N(tt));
    y = linspace(5, 6, N(tt) + 1);
    z = linspace(5, 6, N(tt) - 1);
    
    %% Graded Mesh (dumb octree) 
    weight = 5;
    % Gmesh builds a fat parabola with minimum dictated by wieght
    gmesh = @(x) ((x - (x(end) - x(1))/2 - x(1)).^4 + (x(end) - x(1)) / (length(x) / weight));
    % getdx turns gmesh into a dx and getx does same for x
    getdx = @(x) (x(end) - x(1)) * gmesh(x) / sum(gmesh(x));
    getx = @(x) cumsum([x(1), getdx(x)]);
    % Project cell walls to centres
    centre = @(x) x(1:end - 1) + diff(x);
  
    dz = getdx(z);
    z = centre(getx(z));
    dx = getdx(x);
    x = centre(getx(x));
    dy = getdx(y);
    y = centre(getx(y));

    nz = length(z);
    nx = length(x);
    ny = length(y);
    nnn = nz * nx * ny;
    
    Nnn(tt) = nnn;
      
    [ObsX, ObsY] = meshgrid(5 : 1/10 : 6, 5 : 1/10 : 6);
    ObsZ = zeros(size(ObsX)) + 0.5;
    nobs = numel(ObsX(:));
    
    % Build the dz * dx * dy cell volume vector
    h = dz' * dx;
    h = h(:) * dy;
    h = h(:)';
    
    %% Loop Code
    tic
    g = zeros(nz, nx, ny);
    G = zeros(nobs, nnn);
    count = 1;
    for kk = 1:nobs
        for iy = 1 : ny
            for ix = 1 : nx
                for iz = 1 : nz
                    g(iz, ix, iy) = ((z(iz) - ObsZ(kk) ) * h( (iy - 1) * nz * nx + (ix - 1) * nz + iz)) / ...   % 
                                    ((  (x(ix)-ObsX(kk))^2 +...
                                        (y(iy) - ObsY(kk))^2 +...
                                        (z(iz) - ObsZ(kk))^2 ).^(3/2) );
                end
            end
        end
        G(count, :) = g(:)';
        count = count + 1;
    end
    loopt(tt) = toc; %#ok<*SAGROW>
    
    %% Vectorized Code
    tic
    e = @ (n) ones(1, n);
    % Scale up dimensions to 3D
    kron3 = @(a, b, c)  kron(a, kron(b, c));
    % Kron observations to size G
    kronobs = @(obs) kron( e(nnn), obs(:) );
    % Kron dimensions to size G and subtract obs distance
    H = kron( e(nobs)', h);
    Z = kron( e(nobs)', kron3( e(ny), e(nx), z)) - kronobs(ObsZ);
    X = kron( e(nobs)', kron3( e(ny), x, e(nz))) - kronobs(ObsX);
    Y = kron( e(nobs)', kron3( y, e(nx), e(nz))) - kronobs(ObsY);
    R = (X.^2 + Y.^2 + Z.^2).^(3/2);
    G2 = Z .* H .* 1./R;
    
    vect(tt) = toc;
     
end

fprintf('Norm of (Loop code - Vec code) is %f\n', norm(G - G2))
%figure()
plot(Nnn, loopt, Nnn, vect)
legend('Loop time', 'Vectorized time')
title('Computation Time Comparison of \nVectorized vs Loop based Forward Operator Code')
ylabel('Time [seconds]')
xlabel('Number of cells')