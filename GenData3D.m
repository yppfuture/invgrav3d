 %getpath
%Compute N number of randomely generated ray paths.
%First determine if ray intersects a receiver
%If NO, ray path is only plotted.
%If YES, ray path populate a data matrix, then plotted

close all
clear all
%length of land surveyed
nX=15;
nY=15;
nZ=8;

dz=1;
dx=1;
dy=1;

X0=0;
Y0=0;
Z0=0;

%Define target size [X Y lenght(x) length(y)]
target=[10 10 4];


%Densitiy 2D model
background=0;
anomaly=1;

figure (1)
% load smile;
% load Duo;
% load Trio
% load BigU
% load UpDown
[model]=GenModel3D(nX,nY,nZ,background,anomaly,target);
% model=zeros(nZ,nX,nY);
% model((target(3)-2):(target(3)+2),(target(1)-2):(target(1)+2),(target(2)-2):(target(2)+2))=1;
% imagesc (rho);

hold on

m=reshape(model,nX*nY*nZ,1);


%Create data points
[ObsX,ObsY]=meshgrid(5:12,5:12);
ObsZ=1;


n=size(ObsX,2)*size(ObsY,1);

count=1;
% Build sensitivity matrix
for ii=1:n;
                      
        G(count,:)=GenSenGrav3D(nX,nY,nZ,X0,Y0,Z0,dx,dy,dz,ObsX(ii),ObsY(ii),ObsZ);
count=count+1;
end

set(gca,'YDir','reverse')

%Create data matrix 
data=G*m*1e+6;

%Corrupt with 5% random noise
% d = awgn(data,-12.5);
noise = ((data.*.03).*randn(length(data),1));
d = data + noise;
% d=data.*(unifrnd(-0.05,0.05,length(data),1)+1);

% save ('original','data');
save ('data','data');
save ('kernel','G');
save ('model','m');

% save ('kernel2','G2');

%Map data as a function of receiver and angle

% Plotting original vs corrupted data
d_obs=reshape(d,size(ObsX,2),size(ObsY,1));

figure (1)
imagesc(d_obs)
xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')
 

%Create UBC mesh file without padding cells
save('model.dat','-ascii','m')

fid=fopen('UBC_mesh.msh','w');
fprintf(fid,'%i %i %i\n',nX,nY,nZ)
fprintf(fid,'%i %i %i\n',X0,Y0,Z0)


for jj=1:nY
    fprintf(fid,'%4.2f ',dx);
end
fprintf(fid,'\n',dx)
for ii=1:nX
           fprintf(fid,'%4.2f ',dy);
end
fprintf(fid,'\n',dx)
for kk=1:nZ
       fprintf(fid,'%4.2f ',dz);
end        
fclose(fid);

%Create UBC observation file
count=1;
fid=fopen('UBC_obs.obs','w');
fprintf(fid,'%i\n',length(d));
for ii=1:n

    fprintf(fid,'%4.2f %4.2f %4.2f %e\n',ObsX(ii),ObsY(ii),ObsZ,d(ii)); 
    count=count+1;

end
fclose(fid);