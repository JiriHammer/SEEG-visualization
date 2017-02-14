function brain = getBrainData(plotInfo)
% loads 3D brain MRI or CT = brain
% interpolates
% finds selected channels to plot
% MUST inlcude SPM package in Matlab path

% (c) Jiri, Jan17

%% is SPM12 installed ?
spm_dir = what('spm12');
if isempty(spm_dir) 
    error('SPM12 toolbox required. Install SPM12 and add it to pathdef.m !'); 
end

%% load normalized brain MRI (wT1.nii)
fileName = ['.' filesep plotInfo.file2load];
assert(exist(fileName,'file') == 2);
brain.hdr = spm_vol(fileName);
[brain.vol, brain.xyz] = spm_read_vols(brain.hdr);

%% MNI axis (in mm)
voxSize_old = [nan nan nan];
tmp = reshape(brain.xyz(1,:),size(brain.vol));
x = -squeeze(tmp(:,1,1));
assert(x(1) < x(end));          % ascending order
voxSize_old(1) = x(2) - x(1);

tmp = reshape(brain.xyz(2,:),size(brain.vol));
y = squeeze(tmp(1,:,1))';
assert(y(1) < y(end));
voxSize_old(2) = y(2) - y(1);

tmp = reshape(brain.xyz(3,:),size(brain.vol));
z = squeeze(tmp(1,1,:));
assert(z(1) < z(end));
voxSize_old(3) = z(2) - z(1);

%% interpolation (if needed)
voxSize_new = plotInfo.size_interpolate;     % in [mm]
if voxSize_new == voxSize_old(1) && voxSize_new == voxSize_old(2) && voxSize_new == voxSize_old(3)
    VI = brain.vol;
    xi = x';
    yi = y';
    zi = z';
else
    [X,Y,Z] = meshgrid(x,y,z);
    xi = x(1):voxSize_new:x(end);
    yi = y(1):voxSize_new:y(end);
    zi = z(1):voxSize_new:z(end);
    [XI,YI,ZI] = meshgrid(xi,yi,zi);
    V = permute(brain.vol,[2,1,3]);                        % permutation needed by interp3 (for some reason)
    VI = interp3(X,Y,Z,V,XI,YI,ZI,'spline');
    VI = permute(VI,[2,1,3]);                                   % re-arrange back to [x,y,z] dimensions
end

%% return interpolated volume and axis
assert(size(VI,1) == size(xi,2));
assert(size(VI,2) == size(yi,2));
assert(size(VI,3) == size(zi,2));
brain.VI = VI;
brain.xi = xi;
brain.yi = yi;
brain.zi = zi;
brain.voxSize_new = voxSize_new;
