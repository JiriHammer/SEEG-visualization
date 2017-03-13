function plot_brain3D(vals, plotInfo)
% plots iEEG channel activations on 3D transparent brain model

% (c) Jiri, Jan17

%% select or set here color maps
clrmap.brain = gray(128);                           % colormap for brain (T1, T2, CT, ...)

if isfield(plotInfo, 'colorMap')
    clrmap.chnls = plotInfo.colorMap;
else
    clrmap.chnls = jet(128);                            % default colormap for values: jet
    % clrmap.chnls = getColorMap('bwr', 128);             % colormap for values: blue - white - red
    % clrmap.chnls = getColorMap('bcwwmr', 128);          % colormap for values: blue - cyan - white - magenta - red
end
clrmap.fig = cat(1, clrmap.brain, clrmap.chnls);    % colormap of the figure
alphaVal = 0.2;                                     % transparency of the brain structure (1 = opaque)


%% isosurface of brain grey matter
if ~isfield(plotInfo, 'fv')
    % pass info from loaded brain (see getBrainData.m)
    VI = plotInfo.brain.VI;      % interpolated volume
    xi = plotInfo.brain.xi;      % interpolated x-axis, in [mm] of MNI coors
    yi = plotInfo.brain.yi;      % interpolated y-axis, in [mm] of MNI coors
    zi = plotInfo.brain.zi;      % interpolated z-axis, in [mm] of MNI coors
    assert(size(plotInfo.chnls,2) == size(vals,1));
    V = linTransform(VI, [min(VI(:)), max(VI(:))], [0, 1]);

    separationThreshold = 0.5;                  % (value 0.5 separates gray matter from the dark background). Other values 0 - 1 may work also fine
    fv = isosurface(V, separationThreshold);    % surface, vertex 
else
    fv = plotInfo.fv;
end

%% figure
% --- figure
f = figure('visible','on');
set(f, 'Position', plotInfo.figurePosition);
set(f, 'Colormap', clrmap.fig);
set(f, 'Renderer','OpenGL');
set(f, 'InvertHardcopy','off');                 % preserves black background
set(f, 'PaperPositionMode','auto');
opengl('hardware');

%% channel color limits
if ~isfield(plotInfo, 'chnl_clims')
    clims = [min(vals(:)),max(vals(:))];
else
    clims = plotInfo.chnl_clims;
end

%% colorbar
clrbar_axes = axes('visible','on', 'position',[0.96 0.2 0.01 0.6]);  % position
clrbar_vals = [clims(1):1/1000*diff(clims):clims(2)]';
clrbar_xAx = 1:10;
clrbar_yAx = 1:size(clrbar_vals,1);
clrbar_inds = cVals2cInds(repmat(clrbar_vals,[1,10]), [clims(1),clims(2)], size(clrmap.brain,1)+[1,size(clrmap.chnls,1)]);
clrbar_hndl = image(clrbar_xAx, clrbar_yAx, ones(size(clrbar_inds)));
set(clrbar_hndl, 'CData', clrbar_inds, 'CDataMapping','direct');
minInd = min(clrbar_yAx); maxInd = max(clrbar_yAx); 
tickVals = round([minInd: (maxInd-minInd)/4 :maxInd]);                                                              % tick vals/labels
tickName = cell(1,length(tickVals));
for tick = 1:length(tickVals)
    tickName{tick} = num2str(clrbar_vals(tickVals(tick)), '%01.1f');
end
% set(clrbar_axes, 'XTick',[], 'YTick',tickVals, 'YTickLabel',tickName);
set(gca, 'YDir','normal');
set(gca, 'XTick',[], 'YTick',tickVals, 'YTickLabel',tickName);
set(clrbar_axes, 'YAxisLocation','right');

%% 3D brain model
model_axes = axes('visible','on', 'position',[0.05 0.05 0.9 0.9]);  % position
hold on;
cData = ones(size(fv.vertices,1),1);
brain_cInds = cVals2cInds(cData, [0,1], [1,size(clrmap.brain,1)]);
p = patch('Faces', [fv.faces(:,2),fv.faces(:,1),fv.faces(:,3)], ...
    'Vertices',[fv.vertices(:,2),fv.vertices(:,1),fv.vertices(:,3)], ...
    'EdgeColor','none', 'CData',brain_cInds, 'CDataMapping','direct', 'FaceColor','interp','FaceAlpha',alphaVal);
%isonormals(x,y,z,V,p);

% plot properties
lighting phong
axis equal
% colormap(gray(256))
view(0,90);
camlight headlight
set(gca, 'color', [0 0 0])
%set(gca, 'color', [1 1 1])
set(gca, 'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0])  
view(90,0)
material dull;
set(gca, 'XTick',[], 'YTick',[], 'ZTick',[]);

%% channels
circle_size = 28;
for ch = 1:size(plotInfo.chnls,2)
    i_clr = cVals2cInds(vals(ch), [clims(1),clims(2)], size(clrmap.brain,1)+[1,size(clrmap.chnls,1)]);
    clr = clrmap.fig(i_clr,:);
    [ix,iy,iz] = mni2vox(plotInfo.chnls(ch).MNI_x, plotInfo.chnls(ch).MNI_y, plotInfo.chnls(ch).MNI_z, xi, yi, zi); % index of MNI coor
    scatter3(ix,iy,iz, circle_size, 'MarkerFaceColor',clr, 'MarkerEdgeColor','none');
end

%% output directory
if ~isfield(plotInfo, 'outDir')
    outDir = [pwd filesep '3D_model'];
else
    outDir = plotInfo.outDir;
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end  

%% figure name
if isfield(plotInfo, 'figName')
    figname = plotInfo.figName;
else
    figname = 'notNamed';
end

%% snapshots from different views
for v = 1:size(plotInfo.model_views,1)
    thisView = plotInfo.model_views(v,:);
    view(thisView(1),thisView(2));
    
    % save snapshot
    set(f, 'InvertHardcopy','off');                 % preserves black background
    thisFigName = [figname '_view' num2str(v)];
    if plotInfo.printResolution == 0
        print(f, '-dpng','-r0', [outDir filesep thisFigName '.png']);
    else
        print(f, '-dpng','-r600', [outDir filesep thisFigName '.png']);
    end    
end


%% animation loop (animated gif)
if plotInfo.doAnimation_gif
    outDir_gif = [outDir filesep 'gif_animations'];
    if ~exist(outDir_gif, 'dir')
        mkdir(outDir_gif);
    end  
    filename = [outDir_gif filesep figname '.gif'];
    az = [-180:10:180];
    for n = 1:size(az,2)
        view(az(n),0);
        drawnow;
        frame = getframe(f);
        %frame = getframe(f, get(f, 'Position'));
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if n == 1;
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
    end
end

%% save
saveas(f, [outDir filesep figname '.fig']);
if plotInfo.printResolution == 0
    print(f, '-dpng','-r0', [outDir filesep figname '.png']);
else
    print(f, '-dpng','-r600', [outDir filesep figname '.png']);
end
close(f);    
display(['Figure: ' figname ' stored in: ' outDir]);

