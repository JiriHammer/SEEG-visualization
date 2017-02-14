function rgb_map = getColorMap(cTag, nRows)
% based on 'cTag', returns nRows x 3 RGB colormap

% (c) Jiri, Mar13

%% how many rows?
if nargin < 2
    nRows = 256;
end

%% select colormap
switch cTag
    case 'bwr'
        N = ceil(nRows/4)-1;
        rgb_bw = cat(2, linspace(0,1,N)', linspace(0,1,N)',ones(N,1));
        rgb_wr = cat(2, ones(N,1), linspace(1,0,N)', linspace(1,0,N)');
        rgb_map = cat(1, rgb_bw, rgb_wr);
        
    case 'bwwr'
        rgb_bw = cat(2, [0:255]', [0:255]', 255*ones(256,1))./255;
        rgb_wr = cat(2, 255*ones(256,1), [255:-1:0]', [255:-1:0]')./255;
        rgb_map = cat(1, nthroot(rgb_bw,3), nthroot(rgb_wr,3));
        
    case 'bwwwr'
        x = [0:0.005:1]';
        mu = 0.5;
        sigma = 0.25;   %0.15;
        y = normcdf(x,mu,sigma);    
        N = size(y,1);
        rgb_bw = cat(2, y, y, ones(N,1));
        rgb_wr = cat(2, ones(N,1), flipdim(y,1), flipdim(y,1));
        rgb_map = cat(1, rgb_bw, rgb_wr);        
        
    case 'cwwwm'
        x = [0:0.005:1]';
        mu = 0.5;
        sigma = 0.25;   %0.15;
        y = normcdf(x,mu,sigma);    
        N = size(y,1);
        rgb_cw = cat(2, y, ones(N,1), ones(N,1));
        rgb_wm = cat(2, ones(N,1), flipdim(y,1), ones(N,1));
        rgb_map = cat(1, rgb_cw, rgb_wm);        
        
    case 'wjet'
        rgb_jet = jet(nRows);
        i_cyan = find(rgb_jet(:,1)==0, 1, 'last');
        N = round(i_cyan/2);
        x = [0:N-1]./N;
        mu = 0.5;
        sigma = 0.25;   %0.15;
        y = normcdf(x,mu,sigma);            
        rgb_wc = cat(2, x', ones(length(y),1), ones(length(y),1));
        %rgb_wc = cat(2, y', ones(length(y),1), ones(length(y),1));
        rgb_map = cat(1, flipdim(rgb_wc,1), rgb_jet(i_cyan+1:end,:));

    case 'cwgyr'
        cPoints{1} = [0 1 0 1 1];
        cPoints{2} = [1 1 1 1 0];
        cPoints{3} = [1 1 0 0 0];
        nInterp = round(nRows/(length(cPoints{1})-1));
        rgb_map = nan(nInterp*(length(cPoints{1})-1),3);
        for c = 1:size(cPoints,2)
            cVec = [];
            for s = 1:length(cPoints{c})-1
                cVec = cat(2, cVec, linspace(cPoints{c}(s), cPoints{c}(s+1), nInterp));
            end
            rgb_map(:,c) = cVec';
        end
        
    case 'bcwwmr'
        cPoints{1} = [0.0 0 1 1 1 1.0];
        cPoints{2} = [0.5 1 1 1 0 0.0];
        cPoints{3} = [1.0 1 1 1 1 0.5];
        nInterp = round(nRows/(length(cPoints{1})-1));
        rgb_map = nan(nInterp*(length(cPoints{1})-1),3);
        for c = 1:size(cPoints,2)
            cVec = [];
            for s = 1:length(cPoints{c})-1
                cVec = cat(2, cVec, linspace(cPoints{c}(s), cPoints{c}(s+1), nInterp));
            end
            rgb_map(:,c) = cVec';
        end
        
    case 'b-w-r'
        cPoints{1} = [0 1 1 1];
        cPoints{2} = [0 1 1 0];
        cPoints{3} = [1 1 1 0];
        nInterp = round(nRows/(length(cPoints{1})-1));
        rgb_map = nan(nInterp*(length(cPoints{1})-1),3);
        for c = 1:size(cPoints,2)
            cVec = [];
            for s = 1:length(cPoints{c})-1
                cVec = cat(2, cVec, linspace(cPoints{c}(s), cPoints{c}(s+1), nInterp));
            end
            rgb_map(:,c) = cVec';
        end
        
    case 'mryc'
        rgb_mr = cat(2, 255*ones(256,1), zeros(256,1), [255:-1:0]')./255;
        rgb_ry = cat(2, 255*ones(256,1), [0:255]', zeros(256,1))./255;
        rgb_yg = cat(2, [255:-1:0]', 255*ones(256,1), [0:255]')./255;
        rgb_map = cat(1, rgb_mr, rgb_ry, rgb_yg);        

    case 'myc'
        rgb_my = cat(2, 255*ones(256,1), [0:255]', [255:-1:0]')./255;
        rgb_yc = cat(2, [255:-1:0]', 255*ones(256,1), [0:255]')./255;
        rgb_map = cat(1, rgb_my, rgb_yc);
        
    case 'ycg'
        rgb_gc = cat(2, zeros(256,1), 255*ones(256,1), [0:255]')./255;
        rgb_cy = cat(2, [0:255]', 255*ones(256,1), [255:-1:0]')./255;
        rgb_map = cat(1, rgb_gc, rgb_cy);   
        
    case 'jet'
        rgb_map = jet(nRows);
    case 'hsv'
        rgb_map = hsv(nRows);
    case 'hot'
        rgb_map = hot(nRows);
    case 'copper'
        rgb_map = copper(nRows);
    case 'pink'
        rgb_map = pink(nRows);        
        

    otherwise
        rgb_map = [];
end

% f = figure;
% set(f, 'Colormap', rgb_map);
% imagesc((rand(1000, 1000)-0.5)*2);
% colorbar;
