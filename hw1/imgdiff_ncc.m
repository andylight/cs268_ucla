function [err, varargout] = imgdiff_ncc(I1, I2, varargin)
%IMGDIFF_NCC Computes the Normalized Cross-Correlation between I1 and I2.
%First estimates the 'best' alignment between I1 and I2 via intensity
%differences, then returns the difference. This function will determine
%which image to use as the 'template'.
%Output
%   err = imgdiff_ncc(I1, I2)
%   [err, F] = imgdiff_ncc(I1, I2)
%       F is the NCC response matrix output by normxcorr2.
i_p = inputParser;
i_p.addRequired('I1', @isnumeric);
i_p.addRequired('I2', @isnumeric);
i_p.addParamValue('interactive', false, @islogical);
i_p.parse(I1, I2, varargin{:});
I1 = i_p.Results.I1;
I2 = i_p.Results.I2;
interactive = i_p.Results.interactive;

if size(I1, 1) < size(I2, 1)
    template = I1;
    A = I2;
elseif size(I2, 1) < size(I1, 1)
    template = I2;
    A = I1;
elseif size(I1, 2) < size(I2, 2)    % Heights are equal, check widths
    template = I1;
    A = I2;
elseif size(I2, 2) <= size(I1, 2)
    template = I2;
    A = I1;
end

if size(A, 1) <= size(template, 1)
    diff = size(template, 1) - size(A, 1);
    template = template(1:size(template, 1)-diff, :);
end
if size(A, 2) <= size(template, 2)
    diff = size(template, 2) - size(A, 2);
    template = template(:, 1:size(template, 2)-diff);
end

F = normxcorr2(template, A);
F = F(size(template, 1):end, size(template, 2):end);
[val, loc] = argfind_2d(F, 'max');
% Amount which template should be shifted
%w_Tpad = min(size(template, 2), size(A, 2));
%h_Tpad = min(size(template, 1), size(A, 1));
%yoff = loc(1) - h_Tpad;  xoff = loc(2) - w_Tpad;
yoff = loc(1); xoff = loc(2);

Itmp = nan(size(A));

i0 = yoff + 1; i1 = i0 + size(template, 1) - 1;
j0 = xoff + 1; j1 = j0 + size(template, 2) - 1;

i1 = min([size(template, 1), i1]);
j1 = min([size(template, 2), j1]);
if (i1-i0+1) > size(template, 1)
    disp('HEIGHT MISMATCH');
    4;
elseif (j1-j0+1) > size(template, 2)
    disp('WIDTH MISMATCH');
    5;
end
if i0 <= 0 || j0 <= 0
    disp('wat');
end
Itmp(i0:i1, j0:j1) = template(1:i1-i0+1, 1:j1-j0+1);

%err = nansum(nansum(abs(Itmp - A).^2));
Idiff = abs(Itmp - A).^2;
err = nanmean(Idiff(:));
if nargout == 2
    varargout{1} = F;
end

if interactive == true
    figure;
    subplot(1, 4, 1);
    imshow(A, []);
    title('Image A');
    subplot(1, 4, 2);
    imshow(template, []);
    title('Template');
    subplot(1, 4, 3);
    imshow(F, []);
    title('NCC(A, template)');
    subplot(1, 4, 4);
    Icomp = A;
    Icomp(i0:i1, j0:j1) = template(1:i1-i0+1, 1:j1-j0+1);
    imshow(Icomp, []);
    title(sprintf('Icomposite (err=%.4f loc=(%d, %d))', err, xoff, yoff));
end
end

