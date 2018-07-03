function [brainsurface] = main_brainPlot(vals_channels, mni_channels,names_channels,brainsurface, plotSetup)
%main_brainPlot master script for plotting SEEG values in slices and in 3D model
% function version by Kamil Vlcek 30.5.2018
%
% this script is meant as a tutorial showing a with 3 simple examples
% install: 
%   - add this script to matlab paths (use: Set Paths)
%   - SPM12 (or SPM8) must be installed (add paths to SPM)
% data requirements: 
%   - subject specific MRI (T1.nii) normalized to MNI space
%   - MNI coordinates of individual channels
%   - values of individual channels
% does the following:
%   (A) loads universal brain template (colin 27) and plots data in slices
%   (B) 3D brain model (semi-transparent) with colored channel values
%   (C) loads subject-specific brain data and plots them in slices

% (c) Jiri Hammer, Jan 2017, bug reports: jirihammer@gmail.com
% parametry: vals_channels, mni_channels, names_channels, figureNamePrefix, figureVisible, FontSize, myColorMap,outputDir
%% set path to this script
dir_toolbox = what('SEEG-visualization');     % folder 'visualization_SEEG' MUST be in pathdef.m
if isempty(dir_toolbox) 
    error('Path to toolbox required. Add it to pathdef.m !'); 
end
%addpath(genpath(dir_toolbox.path)); %kamil - path to be saved manually
%dir_curr = pwd;
%cd(dir_toolbox.path);

%assert(exist('outputDir','var')>0,'nedefinovan vystupni adresar');
%outputDir ='d:\eeg\motol\pacienti\p073 Pech VT6\figures\';

niiFile = 'koregistrace_JH\wT1.nii'; %pouziva se jen pri  'SUBJECT-SPECIFIC SLICES'

assert(isfield(plotSetup,'outputDir'), 'outputDir musi byt definovan');
if ~isfield(plotSetup,'figureNamePrefix')
    plotSetup.figureNamePrefix = 'myTest_';
end

if ~isfield(plotSetup,'figureVisible')
    plotSetup.figureVisible = 'on';
end

if ~isfield(plotSetup,'FontSize')
    plotSetup.FontSize = 8; %defaulnti puvodni velikost popisek kanalu
end

if ~isfield(plotSetup,'myColorMap')
    plotSetup.myColorMap = jet(128); %defaultni color scale
end

if ~isfield(plotSetup,'circle_size')
    plotSetup.circle_size = 56; %defaultni velikost kulicky v mozku
end

%% ######################## USER INTERFACE ##################################
% --- interface structure 'plotInfo'
plotInfo = struct(...           % user interface structure: holds most (but not all!) of the user settings
    'outputDir', plotSetup.outputDir, ...   % output directory, where the plots are saved
    'figureNamePrefix',  plotSetup.figureNamePrefix, ...                  % figure name prefix. Program adds automatic suffix to each figure.
    'figurePosition', [1281 -89 1920 964], ...           % position of a figure on screen. For whole screen, evoke a figure, maximize it and type: get(gcf, 'Position')
    'printResolution', 300, ...                           % choices: 0 (= screen resolution) or >0  (= dpi). Resolution of the figures.                          
    'colorScale', [], ... %[0 0.7], ...                               % for example, [-10 10]. If empty, programs adjusts colorscale to 5 & 95 percentile of the data.
    'colorMap',  plotSetup.myColorMap, ...                           % colormap for channel values
    'MRI_fileDir', dir_toolbox.path, ...                % full path to brain MRI NIFTI (.nii) file. Must be normalized to MNI space!
    'size_interpolate', 1.0, ...                        % in [mm], voxel size to which the brain is interpolated; 1.0 means no interpolation
    'size_coloredCube', 3.0, ...                        % in [mm], "voxel" size of the colored channel values
    ... %'slicePlanes', {{'axial','sagittal','coronal'}},... % planes to plot, options: 'axial','sagittal','coronal'
    'slicePlanes', {{'sagittal','axial','coronal'}},... % planes to plot, options: 'axial','sagittal','coronal'
    ... %'model_views', [ 0 89.999; 90 0; 0 0], ...           % camera view angles at which the snapshot of the 3D model are taken. For some odd reason (camera light?), view(0,90) makes grey background...
    ... %'model_views', [90 0], ...           % camera view angles at which the snapshot of the 3D model are taken. For some odd reason (camera light?), view(0,90) makes grey background...
    'doAnimation_gif', false, ...                         % GIF animation of 3D brain model (takes longer time)
    'savefig',false,...                 %if save FIG figures 
    'savepng',true, ...               %if sace png figures
    'circle_size',plotSetup.circle_size,    ....   %size of the circles in 3D plot - 28
    'outputType','3D BRAIN MODEL', ... % ktere obrazky se maji generovat SUBJECT-SPECIFIC SLICES | COLIN27 BRAIN SLICES | 3D BRAIN MODEL
    'niiFileSubject', niiFile, ... 
    'figureVisible', plotSetup.figureVisible, ...
    'FontSize', plotSetup.FontSize, ...
    'customColors', plotSetup.customColors ...
);

% --- load channels MNI coors (variable 'data_channels' in 'channelsInfo.mat')
%load('channelsInfo.mat', 'data_channels');
assert(exist('mni_channels','var') == 1,'MNI coordinates var mni_channels is missing'); %ziskam z headeru pomoci CHHeader.GetMNI
plotInfo.chnls = mni_channels;             % !!! set here your MNI coordinates as a structure array with the same fields as in this example!

if exist('names_channels','var') && ~isempty(names_channels)
   plotInfo.chnames = names_channels;  % jmena kanalu, pokud existuji
end

% --- channels values to plot, format = [channels x time]
%vals_channels = randn(size(data_channels,2),2);   % !!! set here your channel valus, format = [channels x time], in this example: 2 time points of channel values
%vals_channels = ones(size(mni_channels,2),1);    % vsechny elektrody stejnou barvou, jeden casovy okamzik
assert(exist('vals_channels','var')==1,'channel values var vals_channels is missing');

% --- channels color scale
if isempty(plotInfo.colorScale)
    clims = [prctile(vals_channels(:),5), prctile(vals_channels(:),95)]; 
    %clims = max(abs(clims)); 
    if clims(1) == clims(2), clims = sort ( [0 max(vals_channels)] );  end %pokud je napriklad jen jedna hodnota    
    plotInfo.chnl_clims = clims;
else
    plotInfo.chnl_clims = plotInfo.colorScale;
end

if exist('brainsurface','var') && ~isempty(brainsurface)
    plotInfo.fv = brainsurface;  %already computed isosurface - kamil
end 

switch plotInfo.outputType
    
%% ################## (A) COLIN27 BRAIN SLICES ##########################
% colin27-specific: get brain template for slices
    case 'COLIN27 BRAIN SLICES'
        
plotInfo.plottingStyle = 'slices';
plotInfo.MRI_file = [plotInfo.MRI_fileDir filesep 'wT1_subject.nii'];       % subject specific brain normalized to MNI space
plotInfo.brain = getBrainData(plotInfo);

% colin27-specific: plot time-series of brain slices (1 figure / time point)
plotInfo.outDir = [plotInfo.outputDir filesep 'slices_normalizedBrain'];
for t = 1:size(vals_channels,2)          % go through all time points
    plotInfo.figName = [plotInfo.figureNamePrefix 'colin_time' num2str(t)];
    plot_brainSlices(vals_channels(:,t), plotInfo);
end    


%% ################## (B) 3D BRAIN MODEL ##########################
% colin27-specific (but works on other brain's too! Requires segmented grey matter, normalized to MNI space)
% get brain template for 3D model
    case '3D BRAIN MODEL' 

plotInfo.plottingStyle = '3D_model';
plotInfo.size_interpolate = 1;          % consider larger voxels for faster & easier manipulation
plotInfo.MRI_file = [plotInfo.MRI_fileDir filesep  'wc1T1_colin27.nii'];       % gray matter (segmented from colin27 by SPM12)
[plotInfo.brain, plotInfo.fv] = getBrainData(plotInfo); %loads 3D brain MRI or CT = brain
brainsurface = plotInfo.fv; %save for later use



% colin27-specific: plot brain: 3D model 
plotInfo.outDir = [plotInfo.outputDir filesep '3D_model'];
for t = 1:size(vals_channels,2)          % go thru all time points
    disp(['Time point ' num2str(t)]);
    nameSuffix = '';% puvodne 'colin_time' num2str(t) - kamil - nepotrebuju to na nic, kazdy cas ma u me jiny pocet kanaly
    plotInfo.figName = [plotInfo.figureNamePrefix nameSuffix]; 
    plot_brain3D(vals_channels(:,t), plotInfo);
end



%% ################## (C) SUBJECT-SPECIFIC SLICES ##########################
% subject-specific: get brain info for slices
    case 'SUBJECT-SPECIFIC SLICES'
        
plotInfo.plottingStyle = 'slices';
plotInfo.MRI_file = plotInfo.niiFileSubject;       % subject specific brain normalized to MNI space
plotInfo.brain = getBrainData(plotInfo);

% subject-specific: plot time-series of brain slices (1 figure / time point)
plotInfo.outDir = [plotInfo.outputDir filesep 'slices_subjectBrain'];
for t = 1:size(vals_channels,2)          % go thru all time points
    plotInfo.figName = [plotInfo.figureNamePrefix 'subject_time' num2str(t)];
    plot_brainSlices(vals_channels(:,t), plotInfo);
end    


%% change back to previous directory
end
disp('hotovo');
%display('Done. Thanks for trying this tutorial out!  :o) ');
%cd(dir_curr);
