function [Istitch_all, G_all, varargout] = do_stitch_images(imgsdir, ptsdir, varargin)
%DO_STITCH_IMAGES Given images and feature point locations, compute the
%necessary transformations G that correctly 'stitch' together all images.
%
%This function will output multiple Istitch/G pairings if it thinks that
%it found images that belong to different mosaics. This could either occur
%if IMGSDIR contains multiple image mosaics, or the point correspondence
%algorithm failed.
%
%[Istitch_all, G_all] = do_stitch_images(imgsdir, ptsdir, ...)
%   Stitch images in imgsdir with points ptsdir. Let K be the number of
%detected mosaics.
%   Istitch_all is an (1xK) cell array, such that:
%       Istitch_all{i} = Istitch_i
%   G_all is an (1xK) cell array such that:
%       G_all{i} = G_i
%[Istitch_all, G_all, Iblend_all] = do_stitch_images(imgsdir, ptsdir, ...)
%   Also include blended images Iblend_all. Blended images attempt to avoid
%artifacts in the image via various strategies.
%
%Parameters:
%   float 'T' -- Threshold for determining point correspondences, in [0,1].
%                Lower values is stricter, higher is looser.
%   str 'method' -- Image similarity metric to use when determining point
%                   correspondences. One of:
%       'ssd': Sum of Square Differences. Works best.
%       'ncc': Normalized Cross Correlaction. Doesn't work right (TODO).
%   str 'blend' -- Image blending method to reduce artifacts. One of:
%       'overwrite': No blending.
%       'average': When stitching images I0, I1; output average pixel
%                  intensity value for overlapping regions.
%   logical 'show_stitches' -- Display all Istitch in Istitch_all in
%                              separate figures.
%% Parse inputs
i_p = inputParser;
i_p.addRequired('imgsdir', @ischar);
i_p.addRequired('ptsdir', @ischar);
i_p.addParamValue('T', 0.05, @isnumeric);
i_p.addParamValue('method', 'ssd', @ischar);
i_p.addParamValue('blend', 'overwrite', @ischar);
i_p.addParamValue('show_stitches', true, @islogical);
i_p.parse(imgsdir, ptsdir, varargin{:});
imgsdir = i_p.Results.imgsdir;
ptsdir = i_p.Results.ptsdir;
T = i_p.Results.T;
method = i_p.Results.method;
blend = i_p.Results.blend;
show_stitches = i_p.Results.show_stitches;

%% Load images and data
imgpaths = sort(getAllFiles(imgsdir));
imgs = cell([1, length(imgpaths)]); imgs_color = cell([1, length(imgpaths)]);
pts = cell([1, length(imgpaths)]);

for i=1:length(imgpaths)
    imgpath = imgpaths{i};  datapath = get_data_path(imgpath, ptsdir);
    imgs{i} = imread_gray(imgpath, 'range', [0 1]); pts{i} = load(datapath);
    imgs_color{i} = imread_gray(imgpath, 'range', [0 1], 'no_gray', true);
end

%% Compute point correspondences
storedname = sprintf('graph_data_T_%.2f.mat', T);
if exist(storedname, 'file') == 2
    disp(sprintf('[Loading graph, matches from %s]', storedname));
    loadstruct = load(storedname);
    graph = loadstruct.graph;
    matches = loadstruct.matches;
else
    tic;
    disp('Computing point correspondences...');
    [graph, matches] = make_graph(imgs, pts, 'T', T, 'method', method);
    dur_makegraph = toc;
    disp(sprintf('Finished constructing graph (%.4fs)', dur_makegraph));
    save(storedname, 'graph', 'matches', 'T');
end

comps = connected_components(graph);
if numel(comps) > 1
    disp('WARNING - Multiple components detected (output image will not be connected)');
    disp(sprintf('    Nb. components: %d', numel(comps)));
end

%% Compute connected components of the graph
G_all = {};
for gi=1:length(comps)
    comp = comps{gi};
    rootnode = comp{1};
    history = containers.Map('KeyType', 'int32', 'ValueType', 'logical');
    tic;
    G = stitch_graph(rootnode, graph, matches, pts, history);
    G_all{gi} = G;
    dur_stitchgraph = toc;
    disp(sprintf('Finished stitch_graph (%.4fs)', dur_stitchgraph));
end

for gi=1:length(comps)
    G = G_all{gi};  % {{imgidx, G_i}, ...}
    for i=1:length(G)
        G_cur = G{i};
        imgid = G_cur{1};
        T = G_cur{2};
        [theta, sc, tx, ty] = get_trans_params(T);
        theta = rad2deg(theta);
        disp(sprintf('For image %s, trans. parameters are:', imgpaths{imgid}));
        disp(sprintf('    theta=%.4f scale=%.4f tx=%.4f ty=%.4f', theta, sc, tx, ty));
        disp(T);    
    end
end

%% Finally, compute image mosaic for each connected component
Istitch_all = {}; Iblend_all = {};
extents_all = {};   % extents_all{gi}{imgid} := coordinates of pasted image
                    % within reference frame of gi
for gi=1:length(comps)
    G = G_all{gi};
    [x_origin, y_origin] = get_origin(G);
    [wcanvas, hcanvas] = compute_canvas_size(imgs, G);
    wcanvas = 4000; hcanvas = 4000; 
    canvas = zeros([hcanvas wcanvas 3]);
    canvas_blend = zeros([hcanvas wcanvas 3]);
    mask = ones([hcanvas wcanvas]);
    extents = {};
    for i=1:length(G)
        imgid = G{i}{1};
        T = G{i}{2};
        [theta, sc, tx, ty] = get_trans_params(T);
        I = imgs_color{imgid};
        T_rot_scale = [[T(1,1), T(1,2), 0];
                       [T(2,1), T(2,2), 0];
                       [0 0 1]];
        T_trans = inv(T_rot_scale) * T;
        x_offset = T_trans(1,3);
        y_offset = T_trans(2,3);
        % Apply rotation+scale, then translate
        if abs(theta) > 1e-1    % Don't rotate if not necessary
            I = imwarp(I, affine2d(T_rot_scale));
        end
        wI = size(I, 2); hI = size(I, 1);
        %i0 = int32(ty - y_origin + 1); i1 = int32(i0 + hI - 1);
        %j0 = int32(tx - x_origin + 1); j1 = int32(j0 + wI - 1);
        x0 = int32(x_offset - x_origin + 1);
        y0 = int32(y_offset - y_origin + 1);
        if (x0 < 0) || (y0 < 0)
            disp('hi');
        end
        canvas = imgpaste(canvas, I, x0, y0, 'method', 'overwrite');
        [canvas_blend, pt_ul, pt_lr] = imgpaste(canvas_blend, I, x0, y0, 'method', blend, 'mask', mask);
        mask(pt_ul(2):pt_lr(2), pt_ul(1):pt_lr(1)) = 0;
        extents{i} = {pt_ul, pt_lr};
    end
    Istitch_all{gi} = canvas;
    Iblend_all{gi}  = canvas_blend;
    extents_all{gi} = extents;
    if show_stitches
        figure;
        subplot(1, 2, 1); imshow(prepimage(canvas), []);
        title('Raw Image stitch');
        subplot(1, 2, 2); imshow(prepimage(canvas_blend), []);
        title(sprintf('Blended Image Stitch (method=%s)', blend));
        suptitle(sprintf('Image stitch for component=%d/%d\nimgids=%s', gi, length(comps), str_vec(comp)));
    end
end
if nargout == 3
    varargout{1} = Iblend_all;
end
if nargout == 4
    varargout{1} = Iblend_all;
    varargout{2} = extents_all;
end
end

function [theta, sc, tx, ty] = get_trans_params(G)
%GET_TRANS_PARAMS Returns the parameters of this transformation matrix.
theta = atan(G(2,1)/G(1,1));
sc = G(1,1) / cos(theta);

G_rot_scale = [[G(1,1), G(1,2), 0];
               [G(2,1), G(2,2), 0];
               [0 0 1]];
G_trans = inv(G_rot_scale) * G;
tx = G_trans(1,3);
ty = G_trans(2,3);
end

function [x0, y0] = get_origin(Gs)
x0 = 0; y0 = 0;
for i=1:length(Gs)
    G = Gs{i};
    G_mat = G{2};
    [theta, sc, tx, ty] = get_trans_params(G_mat);
    G_rot_scale = [[G_mat(1,1), G_mat(1,2), 0];
                   [G_mat(2,1), G_mat(2,2), 0];
                   [0 0 1]];
    offset = G_rot_scale*[1 1 1]';
    x0_off = offset(1) - 1;
    y0_off = offset(2) - 1;
    x0 = min(x0, abs(sc)*tx + x0_off);
    y0 = min(y0, abs(sc)*ty + y0_off);
end
end

function [w, h] = compute_canvas_size(imgs, G)
%COMPUTE_CANVAS_SIZE Computes the total size spanned by the images in IMGS
%after being transformed by G.
x0 = 0; y0 = 0;
x1 = 0; y1 = 0;
for imgid=1:length(imgs)
    img_i = imgs{imgid};
    T_i = get_T(G, imgid);
    if isnan(T_i)   % imgid is NOT within this component
        continue;
    end
    wI = size(img_i, 2); hI = size(img_i, 1);
    pt0 = T_i * [1 1 1]';   % Upperright corner after T
    pt1 = T_i * [wI hI 1]'; % Lowerright corner after T
    x0 = min([pt0(1), x0]);
    y0 = min([pt0(2), y0]);
    x1 = max([pt1(1), x1]);
    y1 = max([pt1(2), y1]);
end
w = uint32(x1 - x0 + 1);
h = uint32(y1 - y0 + 1);
end

function T = get_T(G, imgid)
%GET_T Returns the transformation matrix T_imgid for a given imgid.
%Unfortunately, entries of G are NOT sorted by imgid (since a particular
%image i may NOT be present in G).
T = nan;
for i=1:length(G)
    g_i = G{i};
    imgid_i = g_i{1};
    T_i = g_i{2};
    if imgid == imgid_i
        T = T_i;
    end
end
end
