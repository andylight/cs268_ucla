function [Imosaic, G] = hw1()
%HW1 Performs the image stitching, outputing the final image mosaic, as
%well as the transformations G_i for each image img_i.
% Imosaic is an RGB image.
% G is a cell array of length N (where N is the number of images) s.t.:
%   G{i} := a cell array of the form:
%       {int imgid, matrix T}
%     where T is the 3x3 transformation matrix that maps image IMGID to
%     the common reference system.
imgsdir = 'imgs';
ptsdir  = 'ptsdata';

% T: Threshold parameter between [0.0, 1.0]. Lower threshold is stricter,
%    higher is looser. Used to determine point correspondences.
T = 0.07;
% METHOD: Which image similarity metric to use.
method = 'ssd';         % Sum-of-Square-Differences
% BLEND: Which image blending strategy to use.
blend  = 'smartcopy';

if exist('getAllFiles') == 0
    addpath_ek();   % Add functions from: ek_util
end

[Istitch_all, G_all, Iblend_all] = do_stitch_images(imgsdir, ptsdir, ...
                                                    'T', T, ...
                                                    'method', method, ...
                                                    'blend', blend);
Istitch = Istitch_all{1};       % Holds the 'naive' image mosaic.
G       = G_all{1};             % Holds the transformation matrices that
                                % transform all images to a common ref.
                                % frame.
Imosaic  = Iblend_all{1};       % Holds the blended (processed) image mosaic.

% Sort G by imgid (increasing order)
keys = cell2mat(cellfun(@(c) c(1), G));
[keys_sorted, idxs] = sort(keys);
G = {G{idxs}};

figure;
imshow(prepimage(Imosaic), []);
suptitle(sprintf('Image Mosaic for %d images in: %s', numel(G), imgsdir));
end
