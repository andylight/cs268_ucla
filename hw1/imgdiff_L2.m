function [err, varargout] = imgdiff_L2(I1, I2, varargin)
%IMGDIFF_L2 Computes the L2 diff (sum-of-square-diffs) between images I1
%and I2. I1 and I2 must be the same dimensions. NaN's will be ignored
%correctly.
Idiff = abs(I1 - I2);
nb_ignored = sum(sum(isnan(Idiff)));
err = nansum(nansum(Idiff)) / (size(Idiff,1)*size(Idiff,2) - nb_ignored);
end

