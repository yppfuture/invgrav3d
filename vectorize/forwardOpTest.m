clear all
close all
% Set up dimensions

N = 6:2:40;
N = 10;

%% Graded Mesh (dumb octree)
curve = 4;
w = 0.1; % Percent finer mesh at finest scale to coarsest.

% Project cell walls to centres
centre = @(x) x(1:end - 1) + 0.5*diff(x);
% Get function for creating meshes with given arguments. See meshfunc for
% details.
getdims = meshfunc(4, 0.1, false);

% Outer Test loop, for testing computational time.
for tt = 1:length(N)
    
    [z, dz] = getdims(linspace(5, 6, N(tt) - 1));
    [x, dx] = getdims(linspace(5, 6, N(tt)));
    [y, dy] = getdims(linspace(5, 6, N(tt) + 1));
    
    zc = centre(z);
    xc = centre(x);
    yc = centre(y);

    nz = length(zc);
    nx = length(xc);
    ny = length(yc);
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
                    g(iz, ix, iy) = ((zc(iz) - ObsZ(kk) ) * h( (iy - 1) * nz * nx + (ix - 1) * nz + iz)) / ...   % 
                                    ((  (xc(ix)-ObsX(kk))^2 +...
                                        (yc(iy) - ObsY(kk))^2 +...
                                        (zc(iz) - ObsZ(kk))^2 ).^(3/2) );
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
    Z = kron( e(nobs)', kron3( e(ny), e(nx), zc)) - kronobs(ObsZ);
    X = kron( e(nobs)', kron3( e(ny), xc, e(nz))) - kronobs(ObsX);
    Y = kron( e(nobs)', kron3( yc, e(nx), e(nz))) - kronobs(ObsY);
    R = (X.^2 + Y.^2 + Z.^2).^(3/2);
    G2 = Z .* H .* 1./R;
    
    vect(tt) = toc;
     
end

fprintf('Norm of (Loop code - Vec code) is %f\n', norm(G - G2))
%figure()
%plot(Nnn, loopt, Nnn, vect)
%legend('Loop time', 'Vectorized time')
%title('Computation Time Comparison of \nVectorized vs Loop based Forward Operator Code')
%ylabel('Time [seconds]')
%xlabel('Number of cells')