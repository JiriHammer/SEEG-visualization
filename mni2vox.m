function [vx, vy, vz] = mni2vox(mx, my, mz, xVector, yVector, zVector)
% converts MNI coors in [mm] to voxel indices

% (c) Jiri, Oct15

vx = closestval(xVector, mx);
vy = closestval(yVector, my);
vz = closestval(zVector, mz);
