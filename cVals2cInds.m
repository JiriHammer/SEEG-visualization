function inds = cVals2cInds(vals, clims, mapInds)
% enables multiple colormaps in 1 figure
% lin. scales 'vals' to 'clims' (values out of range in clims are mapped to their boundaries)
% lin. scales this result to mapInds
% returns indices of the colormap

% (c) Jiri, Feb12, Jan17

assert(clims(1) <= clims(2));
assert(mapInds(1) <= mapInds(2));

if clims(1) == clims(2) && numel(vals)==1 && vals(1) == clims(1)
    inds = 144; %kamil 7.2.2018 - pro pripad, kdy maji vsechny kanaly stejnou hodnoty , napr 0
    %144 je modra barva v clrmap.fig
else
%% truncate values 'vals' to range of 'clims' 

tmp = vals(:);
tmp(tmp<clims(1)) = clims(1);
tmp(tmp>clims(2)) = clims(2);
tmp = reshape(tmp, size(vals));

%% lin. translates this result to mapInds
inds = round(linTransform(tmp, clims, mapInds));

end