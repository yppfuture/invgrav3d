function [LHS] = CompLHS(m, G, nX, nY, nZ)

mcell = nX * nY * nZ;

dmx = zeros(mcell,1);
for ii = 1 : (mcell)
    if ii <= nZ
        dmx(ii) = m(ii) - m(ii + nZ);
    elseif ii > nZ && ii <= (mcell - nX * nZ)
        dmx(ii) = -m(ii - nZ) + 2 * m(ii) - m(ii + nZ);
    elseif ii > (mcell - nX * nZ) && ii <= (mcell - (nX - 1) *nZ)
        dmx(ii) = -m(ii - nZ) + m(ii);
    else
        dmx(ii) = 0;
    end
end

dmy = zeros( mcell, 1);
for ii = 1 : (mcell)
    if ii <= nZ * nX
        dmy(ii) = m(ii) - m(ii + nZ * nX);
    elseif ii > nZ * nX && ii <= (mcell - nX * nZ)
        dmy(ii)= -m(ii - nZ * nX) + 2 * m(ii) - m(ii + nZ * nX);
    else
        dmy(ii)= -m(ii - nZ * nX) + m(ii);
    end
end

dmz = zeros(mcell, 1);
count = 0;
for ii = 1 : mcell
    count = count + 1;
    if count == 1 
       dmz(ii) = m(ii) - m(ii+1);
    elseif count == nZ
        dmz(ii) = -m(ii - 1) + 2 * m(ii);
        count = 0;
    else
        dmz(ii) = -m(ii - 1) + 2 * m(ii) - m(ii + 1);
    end
end

% Use symmetry to speed up.
GtG_rowii = zeros(1,size(G, 2));
GtGm = zeros(mcell, 1);
for ii = 1 : size(G, 2)
    for jj = 1 : size(G,2)
        GtG_rowii(jj) = 0;
        for kk = 1 : size(G, 1)
            GtG_rowii(jj) = GtG_rowii(jj) + G(kk, ii) * G(kk, jj);
        end
    end
    GtGm(ii) = GtG_rowii * m;
end

LHS = dmx + dmy + dmz + GtGm;