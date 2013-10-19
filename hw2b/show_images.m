function show_images(imgsdir, ptsdir)
%SHOW_IMAGES Display images in IMGSDIR, and overlay points from PTSDIR.

imgpaths = sort(getAllFiles(imgsdir));
imgs = cell([1, length(imgpaths)]); imgs_color = cell([1, length(imgpaths)]);
pts = cell([1, length(imgpaths)]);

for i=1:length(imgpaths)
    imgpath = imgpaths{i};  datapath = get_data_path(imgpath, ptsdir);
    imgs{i} = imread_gray(imgpath, 'range', [0 1]); pts{i} = load(datapath);
    imgs_color{i} = imread_gray(imgpath, 'range', [0 1], 'no_gray', true);
end


for i=1:length(imgs)
    I = imgs{i};
    P = pts{i};
    figure;
    subimage(prepimage(I)); hold on;
    for j=1:size(P, 1)
        pt = P(j, :);
        plot(pt(1), pt(2), 'rx'); hold on;
    end
    title(sprintf('Image: %s', imgpaths{i}));
end


end

