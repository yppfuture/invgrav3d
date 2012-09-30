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
[dX,dY,dZ]=getmesh(UBC_mesh);

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
model_comp=conjgrad3D(m0,G,data,nX,nY,nZ);
save('model_out.dat','-ascii','model_comp')
%Model input in compactness
% inv_model_comp=reshape(model_comp,nZ,nX,nY);
figure(1)
imagesc((inv_model_comp))
% title('\bfLeast-Square Solution (for first iteration)')
colormap('jet')
% xlabel('\bfx-axis')
% ylabel('\bfDepth')
caxis([2.7 3.2])
% colorbar
hold on
figure(1);
plot([32 28 28 32 32],[8 8 12 12 8],'r--')    %Single target

% plot([28 25 25 28 28],[8 8 12 12 8],'r--')    %Twins target
% plot([35 32 32 35 35],[8 8 12 12 8],'r--')

% plot([26 24 24 26 26],[8 8 12 12 8],'r-')       %Triplets target
% hold on
% plot([30 28 28 30 30],[3 3 5 5 3],'r-')
% hold on
% plot([37 35 35 37 37],[8 8 11 11 8],'r-')

%Invertion iterations

HS=zeros(1,mcell);
d_error=Wd*(d-G*(m0));
misfit(1)=norm(d_error)^2;
% accuraty(1)=sum(model_comp(rho~=2.7))/sum(model_comp);


WdtWd=Wd'*Wd;


delta=1e-4;    %Damping coefficient
beta=0.0080;         %Trade-off parameter
% beta=[5e-1 5e+1];
lamda=1;       %Weight on compactness
epsilon=1.0;      %Weight on smothness
count=1;
while count<25
%     beta=0.5/(count*10);
%Heaviside step function
    for ii=1:mcell
        HS(ii)=model_comp(ii)/target;
        if HS(ii)<=1.1%||model_comp(i)>.1
            HS(ii)=0;
        else HS(ii)=1;
        end
        Wc(ii,ii)=1/sqrt(model_comp(ii).^2.*(1-HS(ii))+delta);
        
    end
    %get_lamda(G,delta_d,Wc);
    
    delta_d=delta_d-target*G*HS';
    WctWc=(Wc')*Wc;
    
%Weighted data with variance and reference model
% model_comp=(G'*WdtWd*G+beta*(WctWc+20*WtW))\(G'*WdtWd)*delta_d+target*HS';
% beta(rem(count,2)+1)
% Conjugate gradient method
A=(G'*WdtWd*G+beta*Wr'*(lamda*WctWc+epsilon*WtW)*Wr);
RHS=(G'*WdtWd)*delta_d;
model_comp=conjgrad_v3(model_comp-target*HS',A,RHS);

d_obs=G*(model_comp);
d_error=Wd*(d-d_obs);

if norm(d)>norm(d_obs)*10^5
    break
end

count=count+1;
misfit(count)=norm(d_error)^2;
phi_d(count) = norm(Wd*(G*(model_comp)-d))^2;
phi_m(count) =(model_comp)'*(WtW)*(model_comp);


inv_model_comp(:)=model_comp(:)+background;
figure(2)
imagesc((inv_model_comp))
% title('\bfCombined Compact & Smooth')
colormap('jet')
% xlabel('\bfx-axis')
% ylabel('\bfDepth')
caxis([2.7 3.2])
% colorbar

for kk=1:size(Wc,1)
    matWc(kk)=Wc(kk,kk);
end

% figure(4);imagesc(matWc);colorbar

% if misfit(count)>n;
% break
% end

end

% [d_plot]=interp_data(d_obs);



% SNR=-20*log10(norm(d_obs-d)/norm(d));
% figure(6)
% imagesc(receiver,-DIPrange:DIPrange,d_plot)
% % colorbar;
% xlabel('\bfReceiver position')
% ylabel('\bfDip Angle')
% axis equal
% axis tight

figure(2)
hold on
plot([32 28 28 32 32],[8 8 12 12 8],'r--')    %Single target

% plot([24 21 21 24 24],[8 8 12 12 8],'r--')    %Twins target
% plot([39 36 36 39 39],[8 8 12 12 8],'r--')

% plot([26 24 24 26 26],[9 9 12 12 9],'r-')       %Triple target
% hold on
% plot([30 28 28 30 30],[3 3 5 5 3],'r-')
% hold on
% plot([37 35 35 37 37],[8 8 11 11 8],'r-')

figure (3)
[AX,H1,H2] = plotyy(2:count,phi_d(2:end),2:count,phi_m(2:end),'plot');
set(get(AX(1),'Ylabel'),'String','\phi_d') 
set(get(AX(2),'Ylabel'),'String','\phi_m') 
set(H1,'LineStyle','--')
xlabel('Iteration')

figure;plot(d);hold on;plot(d_obs*10^5,'r');



