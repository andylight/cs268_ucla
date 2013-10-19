function out = hw2()
% Computes the transformations between the images w.r.t. the first image
% in the sequence, spirit1983.png.
% Assume a SIMILARITY transformation (translation, rotation, scaling).
%
%   out = hw2()
%       out is a cell matrix of length 3 where:
%           out{i} := {imgpath_i, G_i}
%       where G_i is the transformation warping img_i to spirit1983.png.
imgsdir = 'imgs';
ptsdir  = 'pts';

if exist('getAllFiles') == 0
    addpath_ek();   % Add functions from: ek_util
end

[Istitch_all, G_all, Iblend_all, ~, imgpaths, pts] = do_stitch_images(imgsdir, ptsdir);
G       = G_all{1};             % Holds the transformation matrices that
                                % transform all images to a common ref.
                                % frame.

% Sort G by imgid (increasing order)
keys = cell2mat(cellfun(@(c) c(1), G));
[keys_sorted, idxs] = sort(keys);
G = {G{idxs}};
out = cell([1, length(imgpaths)]);
for i=1:length(G)
    G_cur = G{i};
    imgid = G_cur{1};
    T = G_cur{2};
    [theta, sc, tx, ty] = get_trans_params(T);
    theta = rad2deg(theta);
    disp(sprintf('For image %s, trans. params w.r.t. spirit1983.png are:', imgpaths{imgid}));
    disp(sprintf('    theta=%.4f scale=%.4f tx=%.4f ty=%.4f', theta, sc, tx, ty));
    disp(T);    
    out{i} = {imgpaths{imgid}, T};
end
end
