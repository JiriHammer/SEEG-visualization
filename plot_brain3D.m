function plot_brain3D(vals, plotInfo)
% plots iEEG channel activations on 3D transparent brain model

% (c) Jiri, Jan17

%% select or set here color maps
clrmap.brain = gray(128);                           % colormap for brain (T1, T2, CT, ...)

%% channel color limits
if ~isfield(plotInfo, 'chnl_clims')
    clims = [min(vals(:)),max(vals(:))];
else
    clims = plotInfo.chnl_clims;
end
%%

if plotInfo.customColor
    % pridal 29.6.2018 Kamil
    % custom colormap, pod nulou modra farba, nad nulou cervena
    total = abs(clims(1)-clims(2));
    neg = abs(clims(1));

    blue = ceil((neg/total)*128);
    red = 128 - blue;
    lightblue = ceil(blue/2);
    yellow = ceil(red/2);

    b = [linspace(0.3,1,blue)'; zeros(red,1)];
    g = [zeros(lightblue,1); linspace(0,1,blue-lightblue)'; linspace(1,0,yellow)'; zeros(red-yellow,1)];
    r = [zeros(blue,1) ; linspace(1,0.3,red)'];

    clrmap.chnls = [r g b];

elseif isfield(plotInfo, 'colorMap')
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
    V = linTransform(VI, [min(VI(:)), max(VI(:))], [0, 1]);
    separationThreshold = 0.5;                  % (value 0.5 separates gray matter from the dark background). Other values 0 - 1 may work also fine
    fprintf('computing isosurface for 3D brain model ...' ); %no end of line
    fv = isosurface(V, separationThreshold);    % surface, vertex 
    disp(' done');
else
    fv = plotInfo.fv;
end
xi = plotInfo.brain.xi;      % interpolated x-axis, in [mm] of MNI coors
yi = plotInfo.brain.yi;      % interpolated y-axis, in [mm] of MNI coors
zi = plotInfo.brain.zi;      % interpolated z-axis, in [mm] of MNI coors
assert(size(plotInfo.chnls,2) == size(vals,1),'the number of channels in mni and vals shoudl be the same');

%% figure
% --- figure
f = figure('visible',plotInfo.figureVisible); %kupodivu to funguje, i kdyz je figure invisible - kamil 25.7.2017
if isfield(plotInfo, 'figurePosition')
    set(f, 'Position', plotInfo.figurePosition); %only set position, if defined - kamil
end
set(f, 'Colormap', clrmap.fig);
set(f, 'Renderer','OpenGL');
set(f, 'InvertHardcopy','off');                 % preserves black background
set(f, 'PaperPositionMode','auto');
opengl('hardware');


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


%% output directory
plotInfo.outDir = ['C:\Users\nadabednar\Documents\FGÚ\code\iEEG_scripts' filesep 'slices_subjectBrain'];
if isfield(plotInfo, 'savepng') && plotInfo.savepng==true || isfield(plotInfo, 'savefig') && plotInfo.savefig==true || plotInfo.doAnimation_gif
    %only create dir if needed
    assert(isfield(plotInfo, 'outDir'),'output dir not defined');
    outDir = plotInfo.outDir;    
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end  
end

%% figure name
if isfield(plotInfo, 'figName')
    figname = plotInfo.figName;
else
    figname = 'notNamed';
end

%% snapshots from different views
model_views = zeros(numel(plotInfo.slicePlanes),2); %create view according to named slicePlanes - kamil
for v = 1: numel(plotInfo.slicePlanes)
    if strncmp(plotInfo.slicePlanes{v},'coronal',1) %
        model_views(v,:) = [0 0];
    elseif strncmp(plotInfo.slicePlanes{v},'sagittal',1) %
        model_views(v,:) = [90 0];
    elseif strncmp(plotInfo.slicePlanes{v},'axial',1) %
        model_views(v,:) = [0 89.999]; %For some odd reason (camera light?), view(0,90) makes grey background...
    end
end

for v = 1:size(model_views,1)
    thisView = model_views(v,:);
    view(thisView(1),thisView(2));
    
    % channels - vykreslim je tak, aby byly v popredi - kamil 24.7.2017
    circle_size = plotInfo.circle_size; %28;
    handles_s = zeros(size(plotInfo.chnls,2),2); %druhy sloupec je na names  
    for ch = 1:size(plotInfo.chnls,2)
        i_clr = cVals2cInds(vals(ch), [clims(1),clims(2)], size(clrmap.brain,1)+[1,size(clrmap.chnls,1)]);
        clr = clrmap.fig(i_clr,:);
        [ix,iy,iz] = mnixyz(plotInfo,v,ch,xi,yi,zi,fv,vals(ch),0); %lokalni funkce na souradnice
        handles_s(ch,1) = scatter3(ix,iy,iz, circle_size, 'MarkerFaceColor',clr, 'MarkerEdgeColor','none');
    end
    barvy = {[1 1 0],[1 1 .33],[1 1 .66],[1 1 1]}; %barvy ktere budu stridat u jmen kanalu 
    for ch = 1:size(plotInfo.chnls,2) %jmena kanalu kresilm az potom, aby byla navrchu
        if isfield(plotInfo,'chnames') && (vals(ch)>= mean(clims) || numel(vals<=30) ) %kreslim jen vyrazne kanaly nebo vsechny pokudu jich je malo
            [ix,iy,iz] = mnixyz(plotInfo,v,ch,xi,yi,zi,fv,vals(ch),5); %lokalni funkce na souradnice
            move = rand/2+.5; % mezi 0.5 a 1.0 - aby se napisy neprekryvaly
            if rem(ch,2)==1 %pokud liche cislo %strfind(plotInfo.chnames{ch},'1-') %pokud je jednicka ve jmene kanalu
                %vpravo nad kanalem
                handles_s(ch,2) = text(ix+move,iy+move,iz+move,strrep(plotInfo.chnames{ch},'_','\_'),'FontSize',plotInfo.FontSize,'Color', barvy{rem(ch,4)+1}); %vypisu jmeno kanalu, malym pismem
            else
                %vlevo pod kanalem
                handles_s(ch,2) = text(ix-move,iy-move,iz-move,strrep(plotInfo.chnames{ch},'_','\_'),'FontSize',plotInfo.FontSize,'Color', barvy{rem(ch,4)+1}); %vypisu jmeno kanalu, malym pismem
            end
        end
    end
    % save snapshot
    if isfield(plotInfo, 'savepng') && plotInfo.savepng==true
        set(f, 'InvertHardcopy','off');                 % preserves black background
        thisFigName = [figname '_' plotInfo.slicePlanes{v}]; %kamil - vypisu jmeno pohledu spis nez jeho cislo
        if plotInfo.printResolution == 0
            print(f, '-dpng','-r0', [outDir filesep thisFigName '.png']);
        else
            print(f, '-dpng',['-r' num2str(plotInfo.printResolution)], [outDir filesep thisFigName '.png']);
        end
        display(['Figure: ' thisFigName '.png stored in: ' outDir]);
    end
    for ch = 1:size(plotInfo.chnls,2)
        delete(handles_s(ch,1)); %scatter zase smazu, a nakreslim ho z jineho pohledu
        if handles_s(ch,2) ~= 0
          delete(handles_s(ch,2)); %pokud existuje smazu i text
        end
    end
end


%% animation loop (animated gif)
%tohle ted nebude fungovat, zadne hodnoty kanalu nevykresleny - kamil 24.7.2017
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
        if n == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
    end
end

%% save
%tohle ted nebude fungovat, zadne hodnoty kanalu nevykresleny - kamil 24.7.2017
if isfield(plotInfo, 'savefig') && plotInfo.savefig==true
    saveas(f, [outDir filesep figname '.fig']);
    display(['Figure: ' figname '.fig stored in: ' outDir]);
end
%tohle moc nechapu naco tady je , png uz je ulozeny
if isfield(plotInfo, 'savepng') && plotInfo.savepng==true && size(model_views,1)==0 %model_views already stored, no need to store again
    if plotInfo.printResolution == 0
        print(f, '-dpng','-r0', [outDir filesep figname '.png']);
    else
        print(f, '-dpng',['-r' num2str(plotInfo.printResolution)], [outDir filesep figname '.png']);
    end
    display(['Figure: ' figname '.png stored in: ' outDir]);
end
if isfield(plotInfo, 'savepng') && plotInfo.savepng==true || isfield(plotInfo, 'savefig') && plotInfo.savefig==true
    close(f);    %close figure if it was saved 
end
end
function [ix,iy,iz] = mnixyz(plotInfo,iview,ch,xVector,yVector,zVector,fv,valch,top)
        %kopie z kodu nahore, abych mohl dat nazvy kanalu do popredi pred body
        %kvuli tomu i parametr top, ktery zajisti text pred bodem - 0 pro bod a 5 pro text 
        %fv je isosurface, iview je index v slicePlanes = ktery pohled se vykresluje
        %valch je hodnota eeg tohoto kanalu
        [ix,iy,iz] = mni2vox(plotInfo.chnls(ch).MNI_x, plotInfo.chnls(ch).MNI_y, plotInfo.chnls(ch).MNI_z, xVector, yVector, zVector); % index of MNI coor
        if strncmp(plotInfo.slicePlanes{iview},'coronal',1) %
            iy = min(fv.vertices(:,1))-top-abs(valch); %minimalni y je uplne vepredu u tohoto pohledu
        elseif strncmp(plotInfo.slicePlanes{iview},'sagittal',1) %
            ix = max(fv.vertices(:,2))+top+abs(valch); %maximalni x je uplne vepredu u tohoto pohledu
        elseif strncmp(plotInfo.slicePlanes{iview},'axial',1) %
            iz = max(fv.vertices(:,3))+top+abs(valch); %maximalni z je ulne vepredu u tohoto pohledu
        end
end

