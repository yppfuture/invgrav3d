function [LHS] = compLHS_v3(m, G, Wr, wd, nx, ny, nz)

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
    dmx(ii)=dmx(ii)*Wr(ii);
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
    dmy(ii)=dmy(ii)*Wr(ii);
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
    dmz(ii)=dmz(ii)*Wr(ii);
end

% % Use symmetry to speed up.
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
    LHS(ii) = GtG_mjj + dmx(ii) + dmy(ii) + dmz(ii);
end



return