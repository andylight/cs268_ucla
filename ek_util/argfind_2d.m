function [val, loc] = argfind_2d(mat, fn)
%ARGFIND_2D Implements argmax/argmin.
if strcmp(fn, 'min')
%    [row_vals, row_idxs] = min(mat, [], 1);
%    [col_vals, col_idxs] = min(mat, [], 2);
%    [val, loc_col] = min(row_vals);
%    [val, loc_row] = min(col_vals);
%    loc = [loc_row, loc_col];
    [val, imin] = min(mat(:));
    [i, j] = ind2sub(size(mat), imin(1));
    loc = [i, j];
else
%    [row_vals, row_idxs] = max(mat, [], 1);
%    [col_vals, col_idxs] = max(mat, [], 2);
%    [val, loc_col] = max(row_vals);
%    [val, loc_row] = max(col_vals);
%    loc = [loc_row, loc_col];    
    [val, imax] = max(mat(:));
    [i, j] = ind2sub(size(mat), imax(1));
    loc = [i, j];
end

end
