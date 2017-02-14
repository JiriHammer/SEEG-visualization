%% master script for plotting SEEG values in slices and in 3D model
% (c) Jiri Hammer, Jan 2017, bug reports: jirihammer@gmail.com

%% set paths
dir_toolbox = what('visualization_SEEG');     % folder 'visualization_SEEG' MUST be in pathdef.m
if isempty(dir_toolbox) 
    error('Path to toolbox required. Add it to pathdef.m !'); 
end
addpath(genpath(dir_toolbox.path));
dir_curr = pwd;
cd(dir_toolbox.path);

%% structure 'plotInfo'
plotInfo = struct(...                           % structure that holds most of the data & parameters
    'file2load', 'wT1_subject.nii', ...                 % brain MRI: choices: wT1, wc1T1, wc2T1, ...
    'size_interpolate', 0.5, ...                        % in [mm], voxel size to which the brain is interpolated
    'size_coloredCube', 3.0, ...                        % in [mm], "voxel" size of the colored channel values
    'slicePlanes', {{'axial','sagittal','coronal'}},... % planes to plot, options: 'axial','sagittal','coronal'
    'sliceStep', 1.5 ...                                % in [mm]
    );

%% load channels MNI coors (variable 'data_channels' in 'channelsInfo.mat')
load('channelsInfo.mat', 'data_channels');
assert(exist('data_channels','var') == 1);
plotInfo.chnls = data_channels;

%% channels values to plot
%vals = randn(size(data_channels,2),10);      % format = [channels x time], 10 time points of channel values
vals = randn(size(data_channels,2),2);      % format = [channels x time], 10 time points of channel values

%% channels color limits
clims = [prctile(vals(:),5), prctile(vals(:),95)]; clims = max(abs(clims)); 
plotInfo.chnl_clims = [-clims, clims];

%% subject: get brain info for slices
brain_slices = getBrainData(plotInfo);
plotInfo.brain = brain_slices;

%% subject: plot time-series of brain slices (1 figure / time point)
plotInfo.outDir = [dir_toolbox.path filesep 'figs_subjectsBrain'];
plotInfo.brain = brain_slices;
for t = 1:size(vals,2)          % go thru all time points
    plotInfo.figName = ['subject_' num2str(t)];
    plot_brainSlices(vals(:,t), plotInfo);
end    
            
%% colin27: get brain template for slices
plotInfo.file2load = 'wT1_colin27.nii';
brain_slices = getBrainData(plotInfo);
plotInfo.brain = brain_slices;

%% colin27: plot time-series of brain slices (1 figure / time point)
plotInfo.outDir = [dir_toolbox.path filesep 'figs_normalizedBrain'];
plotInfo.brain = brain_slices;
for t = 1:size(vals,2)          % go thru all time points
    plotInfo.figName = ['colin_' num2str(t)];
    plot_brainSlices(vals(:,t), plotInfo);
end    

%% colin27: get brain template for 3D model
plotInfo.file2load = 'wc1T1_colin27.nii';       % gray matter (segmented from colin27 by SPM12)
params.plot_brainTopo.size_interpolate = 2;     % larger voxels for faster & easier manipulation
brain_model = getBrainData(plotInfo);

%% colin27: plot brain: 3D model 
plotInfo.brain = brain_model;
plotInfo.outDir = [dir_toolbox.path filesep 'figs_normalizedBrain'];
for t = 1:size(vals,2)          % go thru all time points
    plotInfo.figName = ['colin_' num2str(t)];
    plot_brain3D(vals(:,t), plotInfo);
end
