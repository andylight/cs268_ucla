function [Istitch, G] = hw1()
%HW1 Performs the image stitching, outputing the final image mosaic, as
%well as the transformations G_i for each image img_i.
imgsdir = 'imgs';
ptsdir  = 'ptsdata';

T = 0.07;
method = 'ssd';     % Sum-of-Square-Differences
blend  = 'smartcopy';

if exist('getAllFiles') == 0
    addpath_ek();   % Add functions from: ek_util
end

[Istitch_all, G_all, Iblend_all, extents_all] = do_stitch_images(imgsdir, ptsdir, ...
                                                                'T', T, ...
                                                                'method', method, ...
                                                                'blend', blend, ...
                                                                'show_stitches', false);
Istitch = Istitch_all{1};
G       = G_all{1};
Iblend  = Iblend_all{1};
extents = extents_all{1};

% Sort G by imgid (increasing order)
keys = cell2mat(cellfun(@(c) c(1), G));
[keys_sorted, idxs] = sort(keys);
G = {G{idxs}};

figure;
subplot(1, 2, 1); imshow(prepimage(Istitch), []);
title('Raw Image stitch');
subplot(1, 2, 2); imshow(prepimage(Iblend), []); hold on;
for imgid=1:size(G, 2)
    pt_ul = extents{imgid}{1};
    pt_lr = extents{imgid}{2};
    w = pt_lr(1) - pt_ul(1);
    h = pt_lr(2) - pt_ul(2);
    rectangle('Position', [pt_ul(1), pt_ul(2), w, h], 'EdgeColor', 'blue'); hold on;
end
title(sprintf('Blended Image Stitch (method=%s)', blend));
suptitle(sprintf('Image Mosaic for %d images in: %s', numel(G), imgsdir));

imwrite(Iblend, 'mosaic.png');
end
