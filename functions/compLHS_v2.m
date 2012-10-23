function [LHS] = CompLHS(m, G, nx, ny, nz)

mcell = nx * ny * nz;
LHS=zeros( mcell, 1);
dmx = zeros(mcell,1);
for ii = 1 : (mcell)
    if ii <= nz
        dmx(ii) = m(ii) - m(ii + nz);
    elseif ii > nz && ii <= (mcell - nx * nz)
        dmx(ii) = -m(ii - nz) + 2 * m(ii) - m(ii + nz);
    elseif ii > (mcell - nx * nz) && ii <= (mcell - (nx - 1) *nz)
        dmx(ii) = -m(ii - nz) + m(ii);
    else
        dmx(ii) = 0;
    end
end

dmy = zeros( mcell, 1);
for ii = 1 : (mcell)
    if ii <= nz * nx
        dmy(ii) = m(ii) - m(ii + nz * nx);
    elseif ii > nz * nx && ii <= (mcell - nx * nz)
        dmy(ii)= -m(ii - nz * nx) + 2 * m(ii) - m(ii + nz * nx);
    else
        dmy(ii)= -m(ii - nz * nx) + m(ii);
    end
end

dmz = zeros(mcell, 1);
count = 0;
for ii = 1 : mcell
    count = count + 1;
    if count == 1 
       dmz(ii) = m(ii) - m(ii+1);
    elseif count == nz
        dmz(ii) = -m(ii - 1) + 2 * m(ii);
        count = 0;
    else
        dmz(ii) = -m(ii - 1) + 2 * m(ii) - m(ii + 1);
    end
end

% % Use symmetry to speed up.
GtG_rowii = zeros(1,size(G, 2));
% GtGm = zeros(mcell, 1);
for ii = 1 : size(G, 2)
    for jj = ii : size(G, 2)
        GtG_rowii(jj) = 0;
        for kk = 1 : size(G, 1)
            GtG_rowii(jj) = GtG_rowii(jj) + G(kk+size(G,1)*(ii-1)) * G(kk+size(G,1)*(jj-1));
        end
      LHS(ii) = LHS(ii)+GtG_rowii(jj) * m(jj);
      if (ii<size(G, 2))
          LHS(1)=GtG_rowii(1) * m(ii);
      end
    end
    LHS(ii) = LHS(ii)+ dmx(ii) + dmy(ii) + dmz(ii);
end



return