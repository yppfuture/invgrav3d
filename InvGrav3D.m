%Inverting for compact model - Iterative method
%Dominique Fournier
%October 15, 2011

clear all
close all

load data;
load kernel;
load model;
% load original;
UBC_mesh=importdata('UBC_mesh.msh', ' ', 0);
% model=importdata('model.dat', ' ', 0);

%Extract cell dimensions from UBC mesh file
[dX,dY,dZ] = getmesh(UBC_mesh);

nX=length(dX);
nY=length(dY);
nZ=length(dZ);

Lx=min(dX);
Ly=min(dY);
Lz=min(dZ);

% background=2.7; %Background density

mcell=nX*nY*nZ;
mdata=length(data);
% target=0.5;     %Expected density contrast

%Space allocation
% inv_model_QAQC=zeros(size(rho,1),size(rho,2));
% inv_model_comp=zeros(size(rho,1),size(rho,2));

%Modify data and model to substract background density
% model=zeros(size(rho,1)*size(rho,2),1);
% model(:)=rho(:)-background;

m0=zeros(mcell,1);
% m0(:)=background;

%Minimizing the difference between computed and observed data
% delta_data=data;
% delta_d=d;
%% Compute compact model using iterative method
% alphaS=1;
% alphax=alphaS*Lx^2;
% alphaz=alphaS*Lz^2;

alphaS=1/Lx.^2;
alphax=1.0;
alphay=1.0;
alphaz=1.0;

%Build depth weighting
% Wr=get_Wr(nX,nY,nZ);
%Create derivative matrices

% Wx=getWx3D(mcell,nX,nY,nZ);
% Wy=getWy3D(mcell,nX,nY,nZ);
% Wz=getWz3D(mcell,nX,nY,nZ);

% WxtWx=alphax*(Wx')*Wx;
% WztWz=alphaz*(Wz')*Wz;
% dm1x=Wx'*Wx*m;
% dm1y=Wy'*Wy*m;
% dm1z=Wz'*Wz*m;





% dmx=zeros(mcell,1);
% for ii=1:(mcell)
%     if ii<=nZ
%         dmx(ii)=m(ii)-m(ii+nZ);
%     elseif ii>nZ && ii<=(mcell-nX*nZ)
%         dmx(ii)=-m(ii-nZ)+2*m(ii)-m(ii+nZ);
%     elseif ii>(mcell-nX*nZ)&& ii<=(mcell-(nX-1)*nZ)
%         dmx(ii)=-m(ii-nZ)+m(ii);
%     else
%         dmx(ii)=0;
%     end
% end
% 
% dmy=zeros(mcell,1);
% for ii=1:(mcell)
%     if ii<=nZ*nX
%         dmy(ii)=m(ii)-m(ii+nZ*nX);
%     elseif ii>nZ*nX && ii<=(mcell-nX*nZ)
%         dmy(ii)=-m(ii-nZ*nX)+2*m(ii)-m(ii+nZ*nX);
%     else
%         dmy(ii)=-m(ii-nZ*nX)+m(ii);
%     end
% end
% 
% dmz=zeros(mcell,1);
% count=0;
% for ii=1:(mcell)
%     count=count+1;
%     if count==1
%        dmz(ii)=m(ii)-m(ii+1);
%     elseif count==8
%         dmz(ii)=-m(ii-1)+2*m(ii);
%         count=0;
%     else
%         dmz(ii)=-m(ii-1)+2*m(ii)-m(ii+1);
%     end
% end


%Set compactness weighting Wc=I
% Wc=eye(mcell,mcell);
% WctWc=Wc'*Wc;

%Plot compact weighting matrix
% matWc=zeros(m,n);

%Compute data weighting matrix on variance
% Wd=weight_error(data);
% WdtWd=Wd'*Wd;
%Weighted on reference model
% WstWs=alphaS*dx*dz*eye(mcell,mcell);
% WtW=(WxtWx+WztWz+WstWs);
% WtW=WxtWx+WztWz+WstWs;



% Conjugate gradient method for first least square
% A=G'*WdtWd*G+0.075*Wr'*WtW*Wr;
% RHS=G'*WdtWd*delta_d;
model_comp = conjgrad3D(m0,G,data,nX,nY,nZ);
save('model_out.dat','-ascii','model_comp')
