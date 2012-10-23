function [Wr]=CompWr(nx,ny,nz,dx,dy,dz,X0,Y0,Z0,ObsX,ObsY,ObsZ,mode,pow);
% Computes the depth weigthing matrix
% Inputs: nx,ny,nz= Number of cells in East,North and Depth
%         

mcell = nx * ny * nz;
R0 = min( [min(dx) min(dy) min(dz)] );

Wr=zeros(1,mcell);

count=1;
for jj=1:ny
    %First compute de location of the center of the cell
    Y = Y0 + jj * dy(jj) - dy(jj) /2;
    %Then the distance between the cell and the observation
    dY = ( Y - ObsY ) ^2;
    for ii=1:nx
        X = X0 + ii * dx(ii) - dx(ii) /2;
        dX = ( X - ObsX ) ^2;
        for kk=1:nz        
            Z = Z0 + kk * dz(kk) - dz(kk) /2;
            dZ = ( Z - ObsZ ) ^2;
            
            %Then compute the distance between the cell and the observation
            R = ( dX + dY + dZ ) ^(1/2);
            Wr(count) = (1 / (R + R0) ^pow) ^2;
            count = count + 1;
        end
    end
end
