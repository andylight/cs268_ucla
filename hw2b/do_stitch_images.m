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
i_p.addParamValue('rootimgpath', -1);        % Path of image to compute transformations to
i_p.addParamValue('T', 0.05, @isnumeric);
i_p.addParamValue('method', 'ssd', @ischar);
i_p.addParamValue('blend', 'overwrite', @ischar);
i_p.addParamValue('show_stitches', false, @islogical);
i_p.addParamValue('usedisk', false, @islogical);    % If True, save/load graph to disk
i_p.parse(imgsdir, ptsdir, varargin{:});
imgsdir = i_p.Results.imgsdir;
ptsdir = i_p.Results.ptsdir;
rootimgpath = i_p.Results.rootimgpath;
T = i_p.Results.T;
method = i_p.Results.method;
blend = i_p.Results.blend;
show_stitches = i_p.Results.show_stitches;
usedisk = i_p.Results.usedisk;

%% Load images and data
imgpaths = sort(getAllFiles(imgsdir));
imgs = cell([1, length(imgpaths)]); imgs_color = cell([1, length(imgpaths)]);
pts = cell([1, length(imgpaths)]);

if rootimgpath == -1
    % Default to: first imgpath in imgpaths
    rootimgpath = imgpaths{1};
end
rootimgpath = fullfile(rootimgpath);
root_imgid = -1;

for i=1:length(imgpaths)
    imgpath = fullfile(imgpaths{i});
    datapath = get_data_path(imgpath, ptsdir);
    imgs{i} = imread_gray(imgpath, 'range', [0 1]); pts{i} = load(datapath);
    imgs_color{i} = imread_gray(imgpath, 'range', [0 1], 'no_gray', true);
    if strcmp(imgpath, rootimgpath)
        root_imgid = i;
    end
end

if root_imgid == -1
    disp(sprintf('ERROR: Could not find rootimgpath %s', rootimgpath));
    return;
end

%% Compute point correspondences
storedname = sprintf('graph_data_T_%.2f.mat', T);
if exist(storedname, 'file') == 2 && usedisk
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
    if usedisk
        save(storedname, 'graph', 'matches', 'T');
    end
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
    sorted_nodes = sort(cell2mat(comp));
    if ismember(root_imgid, sorted_nodes)
        rootnode = root_imgid;
    else
        % Arbitrarily select first
        rootnode = sorted_nodes(1);
    end
    history = containers.Map('KeyType', 'int32', 'ValueType', 'logical');
    tic;
    G = stitch_graph(rootnode, graph, matches, pts, history);
    G_all{gi} = G;
    dur_stitchgraph = toc;
    disp(sprintf('Finished stitch_graph (%.4fs)', dur_stitchgraph));
end

Istitch_all = {}; Iblend_all = {};
extents_all = {};   % extents_all{gi}{imgid} := coordinates of pasted image
                    % within reference frame of gi
if nargout >= 3
    varargout{1} = Iblend_all;
end
if nargout >= 4
    varargout{2} = extents_all;
end
if nargout >= 5
    varargout{3} = imgpaths;
end
if nargout >= 6
    varargout{4} = pts;
end
if nargout >= 7
    varargout{5} = matches;
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
