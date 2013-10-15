function datapath = get_data_path(imgpath, ptsdir)
%GET_DATA_PATH Returns path to .dat file for input imgpath
[rootdir, fname, ext] = fileparts(imgpath);
datapath = fullfile(ptsdir, sprintf('%s.dat', fname));
if ~exist(datapath)
    datapath = fullfile(ptsdir, sprintf('%s.pts', fname));
end
end
