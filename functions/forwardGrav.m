function [G,Wr] = forwardGrav(nX, nY, nZ, X0, Y0, Z0, dx, dy, dz, ObsX,ObsY,ObsZ)
%Create forward operator

R0 = min( [min(dx) min(dy) min(dz)] );
mcell = nX * nY * nZ;

%Pre-allocate for weighting matrix
Wr=zeros(1,mcell);

%Pre-allocate space for one row of forward operator
G=zeros(1,nX * nY * nZ);

count = 1;

for jj = 1 : nY 
    %First compute de location of the center of the cell
    Y = Y0 + jj * dy(jj) - dy(jj) /2;
    %Then the distance between the cell and the observation
    dY = ( Y - ObsY ) ^2;
   for ii = 1 : nX 
        X = X0 + ii * dx(ii) - dx(ii) /2;
        dX = ( X - ObsX ) ^2;       
       for kk = 1: nZ 
            Z = Z0 - kk * dz(kk) + dz(kk) /2;
            dZ = ( Z - ObsZ ) ^2;
            
            % Compute the distance between the cell and the observation
            R = ( dX + dY + dZ ) ^(1/2);
            Wr(count) = (1 / (R + R0) ^2) ^2;
            
            % Compute the forward operator
            G(count) = (ObsZ - (Z0 - (kk-1) * dz(kk))) * 6.6738e-5 /R ^ 3;
        count = count + 1;
       end
   end
end

clear model;
end
