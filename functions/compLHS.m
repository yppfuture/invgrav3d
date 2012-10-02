function [LHS] = CompLHS(m, G, nx, ny, nz)

mcell = nx * ny * nz;

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

N = size(G, 2);
M = size(G, 1);
GtG_rowii = zeros(1, N);
GtGm = zeros(N, 1);

for ii = 1 : N % Main Loop
    
    for jj = ii : N  % Note loop runs from ii -> N
        GtG_rowii(jj) = 0; % Clear previous results
        
        for kk = 1 : M % Dot product of G columns
            GtG_rowii(jj) = GtG_rowii(jj) + G(kk, ii) * G(kk, jj);
        end     
        GtGm(ii) = GtGm(ii) + GtG_rowii(jj) * m(jj); % Compute top half of triangle
    end
    
    if ii + 1 <= N
        for ll = ii + 1 : N 
            GtGm(ll) = GtGm(ll) + GtG_rowii(ll) * m(ii); % Compute lower half of triangle
        end
    end
end

LHS = dmx + dmy + dmz + GtGm;

return