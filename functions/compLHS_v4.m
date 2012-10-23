function [LHS] = compLHS_v4(m, G, Wr, wd, nX, nY, nZ, dX, dY, dZ, beta)

mcell = nX * nY * nZ;
LHS=zeros( mcell, 1);
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
            dmx(count)=dmx(count)*Wr(count);
            count=count+1;
        end
    end
end
% for ii = 1 : (mcell)
%     if ii <= nZ
%         dmx(ii) = m(ii) - m(ii + nZ);
%     elseif ii > nZ && ii <= (mcell - nX * nZ)
%         dmx(ii) = -m(ii - nZ) + 2 * m(ii) - m(ii + nZ);
%     elseif ii > (mcell - nX * nZ) && ii <= (mcell - (nX - 1) *nZ)
%         dmx(ii) = -m(ii - nZ) + m(ii);
%     else
%         dmx(ii) = 0;
%     end
%     dmx(ii)=dmx(ii)*Wr(ii);
% end

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
            dmy(count)=dmy(count)*Wr(count);
            count=count+1;
        end
    end
end

% for ii = 1 : (mcell)
%     if ii <= nZ * nX
%         dmy(ii) = m(ii) - m(ii + nZ * nX);
%     elseif ii > nZ * nX && ii <= (mcell - nX * nZ)
%         dmy(ii)= -m(ii - nZ * nX) + 2 * m(ii) - m(ii + nZ * nX);
%     else
%         dmy(ii)= -m(ii - nZ * nX) + m(ii);
%     end
%     dmy(ii)=dmy(ii)*Wr(ii);
% end

dmz = zeros(mcell, 1);
count = 1;
% for ii = 1 : mcell
%     count = count + 1;
%     if count == 1 
%        dmz(ii) = m(ii) - m(ii+1);
%     elseif count == nZ
%         dmz(ii) = -m(ii - 1) + 2 * m(ii);
%         count = 0;
%     else
%         dmz(ii) = -m(ii - 1) + 2 * m(ii) - m(ii + 1);
%     end
%     dmz(ii)=dmz(ii)*Wr(ii);
% end
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
            dmz(count)=dmz(count)*Wr(count);
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
            WstWsm(count) = m(count) * dX(ii)* dY(jj) * dZ(kk) * Wr(count);
            count=count+1;
        end
    end
end

%% Compute the misfit term GtG * m
% Use symmetry to speed up.
GtG_rowii = zeros(1,size(G, 2));
GtG_mjj=0;
M=size(G, 2);
N=size(G, 1);

% GtGm = zeros(mcell, 1);
for ii = 1 : M 
    for jj = ii : M
        GtG_rowii(jj) = LHS(jj);
        for kk = 1 : N
            GtG_rowii(jj) = GtG_rowii(jj) + G(kk + N *(ii-1)) * G(kk+ N *(jj-1));
        end
%         GtG_m=GtG_rowii(jj) * m(jj);
      GtG_mjj = GtG_mjj+GtG_rowii(jj) * m(jj) * wd;
      LHS(jj)=LHS(jj)+GtG_rowii(jj) * m(jj) * wd;
%       if (ii<size(G, 2))
%           LHS(1)=GtG_rowii(1) * m(ii);
%       end
    end
    LHS(ii) = GtG_mjj + beta * ( WstWsm(ii) + dmx(ii) + dmy(ii) + dmz(ii) );
end



return