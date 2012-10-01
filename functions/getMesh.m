function [dX, dY, dZ] = getMesh(UBC_mesh)

mesh = UBC_mesh';

xmax = mesh(1);
ymax = mesh(2);
zmax = mesh(3);

header = 6;
for ii = 1 : xmax
    dX(ii) = mesh(header + ii);
end

header = ceil(ii / 3) + header;

for jj = 1 : ymax
    dY(jj) = mesh(header + jj);
end

header = ceil(jj / 3) + header;

for kk = 1 : zmax
    dZ(kk) = mesh(header + kk);
end

end
