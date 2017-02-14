function inds = cVals2cInds(vals, clims, mapInds)
% enables multiple colormaps in 1 figure
% lin. scales 'vals' to 'clims' (values out of range in clims are mapped to their boundaries)
% lin. scales this result to mapInds
% returns indices of the colormap

% (c) Jiri, Feb12, Jan17

assert(clims(1) <= clims(2));
assert(mapInds(1) <= mapInds(2));

%% truncate values 'vals' to range of 'clims' 
tmp = vals(:);
tmp(tmp<clims(1)) = clims(1);
tmp(tmp>clims(2)) = clims(2);
tmp = reshape(tmp, size(vals));

%% lin. translates this result to mapInds
inds = round(linTransform(tmp, clims, mapInds));