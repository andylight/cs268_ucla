function [Istitch_all, G_all] = do_stitch_images(imgsdir, ptsdir, varargin)
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
%
%Parameters:
%   float 'T' -- Threshold for determining point correspondences, in [0,1].
%                Lower values is stricter, higher is looser.
%   str 'method' -- Image similarity metric to use when determining point
%                   correspondences. One of:
%       'ssd': Sum of Square Differences. Works best.
%       'ncc': Normalized Cross Correlaction. Doesn't work right (TODO).
%   logical 'show_stitches' -- Display all Istitch in Istitch_all in
%                              separate figures.
i_p = inputParser;
i_p.addRequired('imgsdir', @ischar);
i_p.addRequired('ptsdir', @ischar);
i_p.addParamValue('T', 0.05, @isnumeric);
i_p.addParamValue('method', 'ssd', @ischar);
i_p.addParamValue('show_stitches', false, @islogical);
i_p.parse(imgsdir, ptsdir, varargin{:});
imgsdir = i_p.Results.imgsdir;
ptsdir = i_p.Results.ptsdir;
T = i_p.Results.T;
method = i_p.Results.method;
show_stitches = i_p.Results.show_stitches;

imgpaths = sort(getAllFiles(imgsdir));
imgs = cell([1, length(imgpaths)]); imgs_color = cell([1, length(imgpaths)]);
pts = cell([1, length(imgpaths)]);

for i=1:length(imgpaths)
    imgpath = imgpaths{i};  datapath = get_data_path(imgpath, ptsdir);
    imgs{i} = imread_gray(imgpath, 'range', [0 1]); pts{i} = load(datapath);
    imgs_color{i} = double(imread(imgpath));
end

storedname = sprintf('graph_data_T_%.2f.mat', T);
if exist(storedname, 'file') == 2
    disp(sprintf('[Loading graph, matches from %s]', storedname));
    loadstruct = load(storedname);
    graph = loadstruct.graph;
    matches = loadstruct.matches;
else
    tic;
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

Istitch_all = {};   G_all = {};

for gi=1:length(comps)
    comp = comps{gi};
    rootnode = comp{1};
    history = containers.Map('KeyType', 'int32', 'ValueType', 'logical');
    tic;
    G = stitch_graph(rootnode, graph, matches, pts, history);
    G_all{gi} = G;
    dur_stitchgraph = toc;
    disp(sprintf('Finished stitch_graph (%.4fs)', dur_stitchgraph));

    [x_origin, y_origin] = get_origin(G);
    [wcanvas, hcanvas] = compute_canvas_size(imgs, G);
    canvas = zeros([hcanvas wcanvas 3]);
    for i=1:length(G)
        curcell = G{i};
        imgid = curcell{1};
        T = curcell{2};
        [theta, tx, ty] = get_rigid_params(T);
        I = imgs_color{imgid};
        if abs(theta) > 1e-1
            T_rot = [[T(1,1), T(1,2), 0];
                    [T(2,1), T(2,2), 0];
                    [0 0 1]];
            I = imwarp(I, affine2d(T_rot));
        end
        wI = size(I, 2); hI = size(I, 1);
        i0 = ty - y_origin + 1; i1 = i0 + hI - 1;
        j0 = tx - x_origin + 1; j1 = j0 + wI - 1;
        i0 = int32(i0); i1 = int32(i1);
        j0 = int32(j0); j1 = int32(j1);
        canvas(i0:i1, j0:j1, 1) = I(:,:,1);
        canvas(i0:i1, j0:j1, 2) = I(:,:,2);
        canvas(i0:i1, j0:j1, 3) = I(:,:,3);
    end
    Istitch_all{gi} = canvas;
    if show_stitches
        figure; imshow(uint8(canvas), []);
        suptitle(sprintf('Image stitch for component=%d/%d\nimgids=%s', gi, length(comps), str_vec(comp)));
    end
end

end

function [theta, tx, ty] = get_rigid_params(G)
theta = acos(min(1, G(1,1)));   % Avoid taking acos(x) for x >= 1 (returns complex values)
tx = G(1, 3); ty = G(2, 3);
end
function [x0, y0] = get_origin(Gs)
x0 = 0; y0 = 0;
for i=1:length(Gs)
    G = Gs{i};
    [theta, tx, ty] = get_rigid_params(G{2});
    x0 = min(x0, tx);
    y0 = min(y0, ty);
end
end

function [w, h] = compute_canvas_size(imgs, G)
%COMPUTE_CANVAS_SIZE Computes the total size spanned by the images in IMGS
%after being transformed by G.
w = 0; h = 0;
for imgid=1:length(imgs)
    img_i = imgs{imgid};
    T_i = get_T(G, imgid);
    if isnan(T_i)   % imgid is NOT within this component
        continue;
    end
    wI = size(img_i, 2); hI = size(img_i, 1);
    pt = T_i * [wI hI 1]'; % Lowerright corner after transformation
    w = uint32(max([pt(1), w]));
    h = uint32(max([pt(2), h]));
end
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
