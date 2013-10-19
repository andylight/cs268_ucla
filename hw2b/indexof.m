function idx = indexof(imgpaths, imgname)
%INDEXOF Find index of imagepath with imgname.
idx = nan;
for i=1:length(imgpaths)
    [~, fname, ~] = fileparts(imgpaths{i});
    if strcmp(fname, imgname)
        idx = i;
        return;
    end
end
end
