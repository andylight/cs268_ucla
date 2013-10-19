function [imgpaths, imgs, ptspaths, pts] = load_data()
%LOAD_DATA Load data.
imgsdir = 'imgs';
ptsdir = 'pts';

imgpaths = sort(getAllFiles(imgsdir));
imgs = cell([1, length(imgpaths)]);
ptspaths = cell([1, length(imgpaths)]);
pts = cell([1, length(imgpaths)]);

for i=1:length(imgpaths)
    imgpath = fullfile(imgpaths{i});
    datapath = get_data_path(imgpath, ptsdir);
    ptspaths{i} = datapath;
    pts{i} = load(datapath);
    imgs{i} = imread_gray(imgpath, 'range', [0 1], 'no_gray', true);
end

end

