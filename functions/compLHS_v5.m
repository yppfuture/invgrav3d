function [LHS] = compLHS_v5(m, G, Wr, wd, nX, nY, nZ, dX, dY, dZ, beta, alphaS)

mcell = nX * nY * nZ;

dmx = zeros(mcell,1);

count=1;

%% Compute the derivative terms ( WxtWx + WytWy + WztWz ) * m
for jj = 1 : nY
    
    for ii = 1 : nX
        
        for kk = 1 : nZ
            if ii == 1
                dmx(count) = m(count) * dX(ii) * dY(jj) * dZ(kk) - ...
                    m(count + nZ) * dX(ii + 1) * dY(jj) * dZ(kk);
                
            elseif ii==nX 
                    dmx(count) = -m(count - nZ) * dX(ii - 1) * dY(jj) * dZ(kk) + ...
                    2*m(count) * dX(ii) * dY(jj) * dZ(kk);
            
            else 
                dmx(count) = -m(count - nZ) * dX(ii - 1)* dY(jj) * dZ(kk) + ...
                    2 * m(count) * dX(ii) * dY(jj) * dZ(kk) - ...
                    m(count + nZ) * dX(ii + 1) * dY(jj) * dZ(kk);
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
                dmy(count) = m(count) * dX(ii) * dY(jj) * dZ(kk) - ...
                    m(count + nZ * nX) * dX(ii) * dY(jj + 1) * dZ(kk);
                
            elseif jj==nY 
                    dmy(count) = -m(count - nZ*nX) * dX(ii) * dY(jj - 1) * dZ(kk) + ...
                    2*m(count) * dX(ii) * dY(jj) * dZ(kk);
            
            else 
                dmy(count) = -m(count - nZ * nX) * dX(ii )* dY(jj - 1) * dZ(kk) + ...
                    2 * m(count) * dX(ii) * dY(jj) * dZ(kk) - ...
                    m(count + nZ* nX) * dX(ii) * dY(jj + 1) * dZ(kk);
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
                dmz(count) = m(count) * dX(ii) * dY(jj) * dZ(kk) - ...
                    m(count + 1) * dX(ii) * dY(jj) * dZ(kk + 1);
                
            elseif kk==nZ 
                    dmz(count) = -m(count - 1) * dX(ii) * dY(jj) * dZ(kk - 1) + ...
                    2*m(count) * dX(ii) * dY(jj) * dZ(kk);
            
            else 
                dmz(count) = -m(count - 1) * dX(ii )* dY(jj) * dZ(kk - 1) + ...
                    2 * m(count) * dX(ii) * dY(jj) * dZ(kk) - ...
                    m(count + 1) * dX(ii) * dY(jj) * dZ(kk + 1);
            end
            dmz(count)=dmz(count) * (1/Wr(count)) ^2;
            count=count+1;
        end
    end
end

%% Compute the smallest term: WstWs * m
% Identity matrix multiplied by length scales is no ref model
WstWsm=zeros(mcell, 1);
count=1;
for jj = 1 : nY
    
    for ii = 1 : nX
        
        for kk = 1 : nZ

            WstWsm(count) = alphaS * m(count) * dX(ii)* dY(jj) * dZ(kk) * (1/Wr(count));

            count=count+1;
        end
    end
end

%% Compute the misfit term GtG * m
% Use symmetry to speed up.
LHS=zeros( mcell, 1);
GtG_mjj=0;
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
    LHS(ii) = LHS(ii) +  ( WstWsm(ii) + dmx(ii) + dmy(ii) + dmz(ii) ) ;
    LHS(ii) = beta * LHS(ii) ;
end



return