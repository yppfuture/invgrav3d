function [model]=GenModel3D(nX,nY,nZ,background,anomaly,target)
%Create density model using Xmax, Zmax,
%density of country rock, density of ore body
%Vector target gives the center position of the ore body
%Size of the target has been fixed to 2 meters across
model=ones(nZ,nX,nY)*background;


for ii=(target(1)+[-2:2])
    for jj=(target(2)+[-2:2])
        for kk=(target(3)+[-2:2])
        model(kk,ii,jj)=anomaly;
        end
    end
end

end