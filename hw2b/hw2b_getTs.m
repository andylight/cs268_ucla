function [G_21, G_32, PtMats] = hw2b_getTs(varargin)
%HW2B_GETTS Compute the (affine) transformations mapping image
%spirit1983.png to spirit1706.png, and the map between spirit1706.png to
%spirit1433.png.
%   [G_21, G_32, PtMats] = hw2b_getTs()
%       G_21 := 3x3 affine matrix mapping spirit1983.png -> spirit1706.png
%       G_32 := 3x3 affine matrix mapping spirit1706.png -> spirit1433.png
%       PtMats := {pts_12, pts_23}, where each pts_ij is a cell array:
%                   {pts_i, pts_j}    pts_i := nx2 matrix of points
%                 such that each point k in pts_i,pts_j correspond
%   spirit1983.png := 1
%   spirit1706.png := 2
%   spirit1433.png := 3
i_p = inputParser;
i_p.addParamValue('K', nan, @isnumeric);
i_p.parse(varargin{:});
K = i_p.Results.K;
imgsdir = 'imgs';
ptsdir  = 'pts_2';

[~, G_all, ~, ~, imgpaths, pts, matches] = do_stitch_images(imgsdir, ptsdir, 'rootimgpath', 'imgs/spirit1983.png', 'K', K);
G       = G_all{1};             % Holds the transformation matrices that
                                % transform all images to a common ref.
                                % frame: spirit1983.png
                                % G{i} := {int imgid, G_imgid}

% Create PtMats
PtMats = create_PtMats(matches, imgpaths, pts);
                                
% Sort G by imgid (increasing order)
keys = cell2mat(cellfun(@(c) c(1), G));
[keys_sorted, idxs] = sort(keys);
G = {G{idxs}};

G_11 = nan;
G_12 = nan;
G_13 = nan;
for i=1:length(G)
    G_cur = G{i};
    imgid = G_cur{1};
    imgpath = imgpaths{imgid};
    T = G_cur{2};    
    [~, filename, ~] = fileparts(imgpath);
    if strcmp(filename, 'spirit1983')
        G_11 = eye(3);
    elseif strcmp(filename, 'spirit1706')
        G_12 = T;
    elseif strcmp(filename, 'spirit1433')
        G_13 = T;
    else
        disp(sprintf('WAT. Unexpected filename: %s', filename));
        return;
    end
end
% Compute G_21, G_32
G_21 = inv(G_12);
G_32 = (inv(G_13) * G_12);
foo = {G_21, G_32};
for i=1:length(foo)
    T = foo{i};
    [theta, sc, tx, ty] = get_trans_params(T);
    theta = rad2deg(theta);
    %disp(sprintf('For image %d, trans. params w.r.t. %d are:', i, i+1));
    %disp(sprintf('    theta=%.4f scale=%.4f tx=%.4f ty=%.4f', theta, sc, tx, ty));
    %disp(T);
end    
end

function PtMats = create_PtMats(matches, imgpaths, pts)
PtMats = cell([1, 2]);
idx1983 = indexof(imgpaths, 'spirit1983');
idx1706 = indexof(imgpaths, 'spirit1706');
idx1433 = indexof(imgpaths, 'spirit1433');
cellthing = matches{idx1983, idx1706};
M1 = cellthing{1};
pts1 = []; pts2 = [];
for i=1:length(M1)
    j = M1(i);
    if ~isnan(j)
        pt1 = pts{idx1983}(i,:);
        pt2 = pts{idx1706}(j,:);
        pts1 = [pts1; pt1];
        pts2 = [pts2; pt2];
    end
end
PtMats{1} = {pts1, pts2};
cellthing = matches{idx1706, idx1433};
M1 = cellthing{1};
pts1 = []; pts2 = [];
for i=1:length(M1)
    j = M1(i);
    if ~isnan(j)
        pt1 = pts{idx1706}(i,:);
        pt2 = pts{idx1433}(j,:);
        pts1 = [pts1; pt1];
        pts2 = [pts2; pt2];
    end
end
PtMats{2} = {pts1, pts2};
end
