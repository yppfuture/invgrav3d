mcell = 100;
nx = 17;
ny = 13;
nz = 11;
dmz = zeros(mcell, 1);
count = 0;
for ii = 1 : mcell
    count = count + 1;
    if count == 1 
        fprintf('%i count == 1 mod = %i\n', ii, mod(ii-1, nz))
    elseif count == nz
        fprintf('%i count == nz mod = %i\n1', ii, mod(ii-1, nz) )
        count = 0;
    else
        fprintf('%i count  else mod = %i\n', ii, mod(ii-1, nz) )
    end
end