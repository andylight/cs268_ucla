function Iout = imgpaste(I, patch, x, y, varargin)
%IMGPASTE Pastes image PATCH onto image I at (X,Y) using various blending
%schemes. (X,Y) refers to the upper-left corner of the destination.
M_OVERWRITE = 'overwrite';
M_AVERAGE   = 'average';
M_MEANSHIFT = 'meanshift';
i_p = inputParser;
i_p.addRequired('I', @isnumeric);
i_p.addRequired('patch', @isnumeric);
i_p.addRequired('x', @isnumeric);
i_p.addRequired('y', @isnumeric);
i_p.addParamValue('method', M_OVERWRITE, @ischar);
i_p.parse(I, patch, x, y, varargin{:});
I = i_p.Results.I;
patch = i_p.Results.patch;
x = i_p.Results.x;
y = i_p.Results.y;
method = i_p.Results.method;

xend = min([size(I, 2), x + size(patch, 2) - 1]);
yend = min([size(I, 1), y + size(patch, 1) - 1]);

is_rgb_I = size(I, 3) > 1;
is_rgb_patch = size(patch, 3) > 1;

% 'Upgrade' images to RGB if necessary
if is_rgb_I && ~is_rgb_patch
    patch = gray2rgb(patch);
elseif ~is_rgb_I && is_rgb_patch
    I = gray2rgb(I);
end
Iout = I;   % Copies matrix I (Matlab uses copy-on-write)
if strcmp(method, M_OVERWRITE)
    if is_rgb_I
        Iout(y:yend, x:xend, :) = patch;
    else
        Iout(y:yend, x:xend) = patch;
    end
elseif strcmp(method, M_AVERAGE)
    if is_rgb_I
        Iout_patch = Iout(y:yend, x:xend, :);
        patch_blended = (Iout_patch + patch) / 2;
        Iout(y:yend, x:xend, :) = patch_blended;
    else
        Iout_patch = Iout(y:yend, x:xend);
        patch_blended = (Iout_patch + patch) / 2;
        Iout(y:yend, x:xend) = patch_blended;
    end
elseif strcmp(method, M_MEANSHIFT)
    % Assume that I, patch differ only by a constant difference, i.e.:
    %   I = patch + t;
    %xstart = nanmin([nanmax([firstnan_row(Iout, x, y), x]), x+size(patch,2)-1]);
    %ystart = nanmin([nanmax([firstnan_col(Iout, x, y), y]), y+size(patch,1)-1]);
    if is_rgb_I
        Ioverlap = Iout(y:nanmin([yend, firstnan_col(Iout, x, y)]), ...
                        x:nanmin([xend, firstnan_row(Iout, x, y)]), ...
                        :);
        mean1_r = mean2(Ioverlap(:,:,1));
        mean1_g = mean2(Ioverlap(:,:,2));
        mean1_b = mean2(Ioverlap(:,:,3));
        mean2_r = mean2(patch(:,:,1));
        mean2_g = mean2(patch(:,:,2));
        mean2_b = mean2(patch(:,:,3));
        patch_shift = patch;
        if ~isnan(mean1_r)
            patch_shift(:,:,1) = patch_shift(:,:,1) + (mean1_r - mean2_r);
            patch_shift(:,:,2) = patch_shift(:,:,2) + (mean1_g - mean2_g);
            patch_shift(:,:,3) = patch_shift(:,:,3) + (mean1_b - mean2_b);        
            patch_shift(patch_shift > 1) = 1;
            patch_shift(patch_shift < 0) = 0;
        end
        Iout(y:yend, x:xend, :) = patch_shift;        
        %Iout(ystart:yend, xstart:xend, :) = patch_shift((ystart-y+1):end, (xstart-x+1):end, :);
        %mask = isnan(Iout(y:yend, x:xend, :));
        
    else
        Ioverlap = Iout(y:min([yend, firstnan_row(Iout, x, y)]), ...
                        x:min([xend, firstnan_col(Iout, x, y)]));
        meanA = mean2(Ioverlap);
        meanB = mean2(patch);
        if ~isnan(meanA)
            patch_shift(:,:) = patch(:,:) + (meanA - meanB);
        end
        Iout(y:yend, x:xend) = patch_shift;
    end    
end
end

function Iout = gray2rgb(I)
%GRAY2RGB Outputs single-channel image I to a three-channel RGB image.
if size(I, 3) > 1
    Iout = I;
    return;
end
Iout =  zeros([size(I), 3]);
Iout(:,:,1) = I;
Iout(:,:,2) = I;
Iout(:,:,3) = I;
end

function out = firstnan_row(I, x, y)
%FIRSTNAN_ROW Find first-nan value row-wise in I, starting at (x,y).
out = nan;
for j=x:size(I, 2)
    if isnan(I(y,j))
        out = j;
        return;
    end
end
end

function out = firstnan_col(I, x, y)
%FIRSTNAN_COL Find first-nan value col-wise in I, starting at (x,y).
out = nan;
for i=y:size(I, 1)
    if isnan(I(i,x))
        out = i;
        return;
    end
end
end
