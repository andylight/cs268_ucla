function patch = get_img_patch(I, x, y, w, h, varargin)
%GET_IMG_PATCH Returns patch from I centered at (x,y) with dims (2*w,2*h). If
%the patch extends outside I, zero-padding is added.
i_p = inputParser;
i_p.addRequired('I', @isnumeric);
i_p.addRequired('x', @isnumeric);   i_p.addRequired('y', @isnumeric);
i_p.addRequired('w', @isnumeric);   i_p.addRequired('h', @isnumeric);
i_p.addParamValue('fillval', 0, @isnumeric);
i_p.parse(I, x, y, w, h, varargin{:});
I = i_p.Results.I; 
x = i_p.Results.x;  y = i_p.Results.y;
w = i_p.Results.w;  h = i_p.Results.h;
fillval = i_p.Results.fillval;

wimg = size(I, 2);    himg = size(I, 1);
patch = zeros([(2*h)+1 (2*w)+1]);
if fillval ~= 0
    patch(:,:) = fillval;
end
y0 = max([1, y - h]);
y1 = min([himg, y + h]);
x0 = max([1, x - w]);
x1 = min([wimg, x + w]);

i0 = abs(min([0, y-h-1])) + 1;
i1 = i0 + (y1-y0);
j0 = abs(min([0, x-w-1])) + 1;
j1 = j0 + (x1-x0);
%patch(1:(y1-y0+1), 1:(x1-x0+1)) = I(y0:y1, x0:x1);
patch(i0:i1, j0:j1) = I(y0:y1, x0:x1);
end
