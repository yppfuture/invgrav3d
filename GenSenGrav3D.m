function [G]=GenSenGrav3D(nX,nY,nZ,X0,Y0,Z0,dx,dy,dz,ObsX,ObsY,ObsZ)
%Create sensitivity matrix

%Pre-allocate space for one row
G=zeros(1,nX*nY*nZ);

count=1;

for jj=0:nY-1
   for ii=0:nX-1
       for kk=0:-1:-(nZ-1)

   G(count)=6.6738e-11/(((ObsX-(X0+ii*dx))^2+(ObsY-(Y0+jj*dy))^2+(ObsZ-(Z0+kk*dz))^2));
   count=count+1;
       end
   end
   
end

clear model;
end