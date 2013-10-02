function datapath = get_data_path(imgpath, ptsdir)
%GET_DATA_PATH Returns path to .dat file for input imgpath
[rootdir, fname, ext] = fileparts(imgpath);
split0 = strsplit(fname, '_');
prefix = split0{1};
idx = split0{2};
datapath = fullfile(ptsdir, sprintf('%s_%s.dat', prefix, idx));
end
