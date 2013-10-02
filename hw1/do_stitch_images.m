function [Istitch, G] = do_stitch_images(imgsdir, ptsdir, varargin)
%DO_STITCH_IMAGES Given images and feature point locations, compute the
%necessary transformations G that correctly 'stitch' together all images.
i_p = inputParser;
i_p.addRequired('imgsdir', @ischar);
i_p.addRequired('ptsdir', @ischar);
i_p.addParamValue('T', 0.05, @isnumeric);
i_p.parse(imgsdir, ptsdir, varargin{:});
imgsdir = i_p.Results.imgsdir;
ptsdir = i_p.Results.ptsdir;
T = i_p.Results.T;

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
    [graph, matches] = make_graph(imgs, pts, 'T', T);
    dur_makegraph = toc;
    disp(sprintf('Finished constructing graph (%.4fs)', dur_makegraph));
    save(storedname, 'graph', 'matches', 'T');
end

comps = connected_components(graph);
if numel(comps) > 1
    disp('WARNING - Multiple components detected (output image will not be connected)');
    disp(sprintf('    Nb. components: %d', numel(comps)));
end

canvases = {};  % Stores all image stitches for each connected component

for gi=1:length(comps)
    comp = comps{gi};
    rootnode = comp{1};
    history = containers.Map('KeyType', 'int32', 'ValueType', 'logical');
    tic;
    G = stitch_graph(rootnode, graph, matches, pts, history);
    dur_stitchgraph = toc;
    disp(sprintf('Finished stitch_graph (%.4fs)', dur_stitchgraph));

    [x_origin, y_origin] = get_origin(G);
    canvas = zeros([1000 1000 3]);
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
    canvases{gi} = canvas;
    figure; imshow(uint8(canvas), []);
    suptitle(sprintf('Image stitch for component=%d/%d\nimgids=%s', gi, length(comps), str_vec(comp)));
end

Istitch = canvases{1}; % Arbitrarily output the first image stitch

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
    T_i = G{imgid};
    wI = size(img_i, 2); hI = size(img_i, 1);
    pt = T_i * [wI hI 1]'; % Lowerright corner after transformation
    w = max([pt(1), w]); h = max([pt(2), h]);
end
end
    