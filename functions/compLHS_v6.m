function [LHS] = compLHS_v6(m, G, Wr, WctWc, wd, nX, nY, nZ, dX, dY, dZ, beta, alphaS)

mcell = nX * nY * nZ;

dmx = zeros(mcell,1);

count=1;

%% Compute the derivative terms ( WxtWx + WytWy + WztWz ) * m
for jj = 1 : nY
    
    for ii = 1 : nX
        
        for kk = 1 : nZ
            if count <= nZ
                dmx(count) = (m(count)  - ...
                    m(count + nZ) ) * dY(jj) * dZ(kk) / dX(ii);
                
            elseif count > nZ &&  count <= ( nY * nZ * (nX-1))
                    dmx(count) = ( -m(count - nZ) + -m(count + nZ) + ...
                    2 * m(count) ) * dY(jj) * dZ(kk) / dX(ii);
            
            elseif  count > ( nY * nZ * (nX-1)) && count <= ( nY * nZ * (nX-1))+ nZ
                dmx(count) = ( m(count)  - ...
                    m(count - nZ)) * dY(jj) * dZ(kk) / dX(ii);
            end
            dmx(count)=dmx(count) * (1/Wr(count));
            count=count+1;
        end
    end
end


dmy = zeros( mcell, 1);
count=1;
for jj = 1 : nY
    
    for ii = 1 : nX
        
        for kk = 1 : nZ
            if jj == 1
                dmy(count) = (m(count) - ...
                    m(count + nZ * nX) ) * dX(ii) * dZ(kk) / dY(jj);
                
            elseif jj==nY 
                    dmy(count) = (-m(count - nZ*nX) + ...
                    m(count)) * dX(ii) * dZ(kk) / dY(jj);
            
            else 
                dmy(count) = (-m(count - nZ * nX) + ...
                    2 * m(count) - ...
                    m(count + nZ* nX) ) * dX(ii) * dZ(kk) / dY(jj);
            end
            dmy(count)=dmy(count) * (1/Wr(count));
            count=count+1;
        end
    end
end


dmz = zeros(mcell, 1);
count = 1;

for jj = 1 : nY
    
    for ii = 1 : nX
        
        for kk = 1 : nZ
            if kk == 1
                dmz(count) = (m(count) - ...
                    m(count + 1)) * dX(ii) * dY(jj) / dZ(kk);
                
            elseif kk==nZ 
                    dmz(count) = (-m(count - 1) + ...
                    m(count)) * dX(ii) * dY(jj) / dZ(kk);
            
            else 
                dmz(count) = (-m(count - 1) + ...
                    2 * m(count) - ...
                    m(count + 1)) * dX(ii) * dY(jj) / dZ(kk);
            end
            dmz(count)=dmz(count) * (1/Wr(count)) ^2;
            count=count+1;
        end
    end
end

%% Compute the smallest term: WstWs * m
% Identity matrix multiplied by length scales is no ref model
WstWsm=zeros(mcell, 1);
WctWcm=zeros(mcell, 1);
count=1;
for jj = 1 : nY
    
    for ii = 1 : nX
        
        for kk = 1 : nZ

            WstWsm(count) = alphaS * m(count) * dX(ii)* dY(jj) * dZ(kk) * (1/Wr(count));
            WctWcm(count) = WctWc (count) * m(count);
            count=count+1;
        end
    end
end

% %% Compute the compact term: WctWc * m
% % Identity matrix multiplied by length scales is no ref model
% WctWcm=zeros(mcell, 1);
% count=1;
% alphaC=10;
% for jj = 1 : nY
%     
%     for ii = 1 : nX
%         
%         for kk = 1 : nZ
% 
%             WctWcm(count) = alphaC * dX(ii)* dY(jj) * dZ(kk) * m(count) / ( m(count) ^2 + 1e-10) ^ (0.55);
% 
%             count=count+1;
%         end
%     end
% end
% 
% WctWcm = WctWcm / max(WctWcm)* max(m);

%% Compute the misfit term GtG * m
% Use symmetry to speed up.
LHS=zeros( mcell, 1);
M=size(G, 2);
N=size(G, 1);

% GtGm = zeros(mcell, 1);
for ii = 1 : M 

    for jj = ii : M        
        GtG_rowii=0;

        for kk = 1 : N
            GtG_rowii = GtG_rowii + G(kk + N *(ii-1)) * G(kk+ N *(jj-1)) * wd(kk) ;
        end

        
      LHS(ii) = LHS(ii) + GtG_rowii * m(jj)  ;

          if (ii~=jj)
              LHS(jj) = LHS(jj) + GtG_rowii * m(ii)   ;
          end


    end
    LHS(ii) = LHS(ii) +  ( dmx(ii) + dmy(ii) + dmz(ii) + WctWcm (ii) ) ;
    LHS(ii) = beta * LHS(ii) ;
end



return