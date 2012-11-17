 %getpath
%Compute N number of randomely generated ray paths.
%First determine if ray intersects a receiver
%If NO, ray path is only plotted.
%If YES, ray path populate a data matrix, then plotted

close all
clear all
addpath data
addpath functions





% For now all the cells have dimension 1x1x1
dX = [1 1 ones(1,12) 1 1];
dY = [1 1 ones(1,12) 1 1];
dZ = [ones(1,10) 1 1];

X0 = 0;
Y0 = 0;
Z0 = 0;

%length of land surveyed
nX = length(dX);
nY = length(dY);
nZ = length(dZ);


mcell=nX*nY*nZ;

%Densitiy 2D model
background = 0;
anomaly = 0.5;

figure (1)
% load smile;
% load Duo;
% load Trio
% load BigU
% load UpDown
target = [10 12 5];
[model] = genModel(nX, nY, nZ, background, anomaly, target);

target = [4 5 6];
[model] = genModel(nX, nY, nZ, model, anomaly, target);

target = [12 5 5];
[model] = genModel(nX, nY, nZ, model, anomaly, target);

% model=zeros(nZ,nX,nY);
% model((target(3)-2):(target(3)+2),(target(1)-2):(target(1)+2),(target(2)-2):(target(2)+2))=1;
% imagesc (rho);

hold on

m=reshape(model,nX*nY*nZ,1);


%Create data points
[ObsX, ObsY] = meshgrid(2:1:14,2:1:14);
ObsZ = 1;
n = size(ObsX, 2) * size(ObsY, 1);

count = 1;

% Initiate G matrix
G=zeros(n,mcell);

% Compute depth weigthing matrix
% mode 0: distance weigthing , 1: depth weighting
% pow: Power of the exponentiel decay (default 2 for grav, 3 for mag)
Wr=zeros(1,mcell);

% Build forward operator
for ii=1:n;
       [G(count,:),Wr] = forwardGrav_v2(nX, nY, nZ, X0, Y0, Z0,...
            dX, dY, dZ, ObsX(ii), ObsY(ii), ObsZ);
%         wr=CompWr(nX,nY,nZ,dX,dY,dZ,X0,Y0,Z0,ObsX(ii), ObsY(ii), ObsZ, mode,pow);
%         Wr=Wr+wr;
%         G(count,:)=G(count,:) .* (1 ./ wr) .^2;
count = count + 1;
end

% Square root for the sum of the squares
% Plus another square root of result because inside the objective function, 
% but we will square after in WrtWr ... so only one sqrt.
% Wr=Wr.^(1/2);

% Normalize depth weighting with the largest value
Wr=Wr./max(Wr);

% % Re-weight the kernel with distance weigthing function
% for ii= 1 : size (G,1)
%     G(ii,:)=G(ii,:) .* (1 ./ (Wr.^1));
% end

%Create data matrix
data = G * m ;

%Corrupt with 5% random noise
% d = awgn(data,-12.5);
noise = ( (data.*.02) .* randn(length(data),1) );
d = data + noise;
% d=data.*(unifrnd(-0.05,0.05,length(data),1)+1);

% save ('original','data');
save('data/data.mat','d');
save('data/kernel.mat','G');
save('data/model.mat','m');
save('data/Wr.mat','Wr');

% save ('kernel2','G2');

d_obs = reshape(d, size(ObsY,1), size(ObsX,2));

figure (1)
imagesc(d_obs)
xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')

%Create UBC mesh file without padding cells
save('data/model.dat', '-ascii', 'm')

%Create UBC observations file 
save('data/Obs_grav.dat', '-ascii', 'd')

%Create UBC model weight file 
Wr_out=Wr';
save('data/Wr.dat', '-ascii', 'Wr_out')

fid=fopen('data/UBC_mesh.msh', 'w');
fprintf(fid, '%i %i %i\n', nX, nY, nZ);
fprintf(fid, '%i %i %i\n', X0, Y0, Z0);


for jj=1:nX
    fprintf(fid, '%4.2f ', dX(jj));
end
fprintf(fid,'\n');
for ii=1:nY
           fprintf(fid,'%4.2f ', dY(ii));
end
fprintf(fid, '\n');
for kk=1 : nZ
       fprintf(fid, '%4.2f ', dZ(kk));
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
