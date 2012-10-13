 %getpath
%Compute N number of randomely generated ray paths.
%First determine if ray intersects a receiver
%If NO, ray path is only plotted.
%If YES, ray path populate a data matrix, then plotted

close all
clear all
addpath data
addpath functions

%length of land surveyed
nx = 15;
ny = 15;
nz = 8;

mcell=nx*ny*nz;

% For now all the cells have dimension 1x1x1
dx = ones(1,nx);
dy = ones(1,ny);
dz = ones(1,nz);

X0 = 0;
Y0 = 0;
Z0 = 0;

%Define target size [X Y lenght(x) length(y)]
target = [10 10 6];


%Densitiy 2D model
background = 0;
anomaly = 500;

figure (1)
% load smile;
% load Duo;
% load Trio
% load BigU
% load UpDown
[model] = genModel(nx, ny, nz, background, anomaly, target);
% model=zeros(nz,nx,ny);
% model((target(3)-2):(target(3)+2),(target(1)-2):(target(1)+2),(target(2)-2):(target(2)+2))=1;
% imagesc (rho);

hold on

m=reshape(model,nx*ny*nz,1);


%Create data points
[ObsX, ObsY] = meshgrid(5:12,5:12);
ObsZ = 1;
n = size(ObsX, 2) * size(ObsY, 1);

count = 1;

% Initiate G matrix
G=zeros(n,mcell);

% Compute depth weigthing matrix
% mode 0: distance weigthing , 1: depth weighting
% pow: Power of the exponentiel decay (default 2 for grav, 3 for mag)
Wr=zeros(1,mcell);
mode=0;
pow=2;

% Build forward operator
for ii=1:n;
       [G(count,:),wr] = forwardGrav(nx, ny, nz, X0, Y0, Z0,...
            dx, dy, dz, ObsX(ii), ObsY(ii), ObsZ);
%         wr=CompWr(nx,ny,nz,dx,dy,dz,X0,Y0,Z0,ObsX(ii), ObsY(ii), ObsZ, mode,pow);
        Wr=Wr+wr;
count = count + 1;
end

% Square root for the sum of the squares
% Plus another square root of result because inside the objective function, 
% but we will square after in WrtWr ... so only one sqrt.
Wr=Wr.^(1/2);

% Normalize depth weighting with the largest value
Wr=Wr./max(Wr);

%Create data matrix
data = G * m ;

%Corrupt with 5% random noise
% d = awgn(data,-12.5);
noise = ( (data.*.03) .* randn(length(data),1) );
d = data + noise;
% d=data.*(unifrnd(-0.05,0.05,length(data),1)+1);

% save ('original','data');
save('data/data.mat','data');
save('data/kernel.mat','G');
save('data/model.mat','m');
save('data/Wr.mat','Wr');

% save ('kernel2','G2');

d_obs = reshape(d, size(ObsX,2), size(ObsY,1));

figure (1)
imagesc(d_obs)
xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')

%Create UBC mesh file without padding cells
save('data/model.dat', '-ascii', 'm')

%Create UBC observations file 
save('data/Obs_grav.dat', '-ascii', 'd')

fid=fopen('data/UBC_mesh.msh', 'w');
fprintf(fid, '%i %i %i\n', nx, ny, nz)
fprintf(fid, '%i %i %i\n', X0, Y0, Z0)


for jj=1:nx
    fprintf(fid, '%4.2f ', dx(jj));
end
fprintf(fid,'\n',dx)
for ii=1:ny
           fprintf(fid,'%4.2f ', dy(ii));
end
fprintf(fid, '\n', dx)
for kk=1 : nz
       fprintf(fid, '%4.2f ', dz(kk));
end
fclose(fid);

%Create UBC observation file
count = 1;
fid = fopen('data/UBC_obs.obs','w');
fprintf(fid,'%i\n',length(d));
for ii=1:n

    fprintf(fid,'%4.2f %4.2f %4.2f %e\n',...
        ObsX(ii), ObsY(ii), ObsZ, d(ii));
    count = count + 1;

end
fclose(fid);
