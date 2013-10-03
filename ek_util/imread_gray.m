function I = imread_gray(imgpath, varargin)
%IMREAD_GRAY Loads in an image as grayscale, double-precision.
i_p = inputParser;
i_p.addRequired('imgpath', @ischar);
i_p.addParamValue('range', '');  % [LOW HIGH] or '' for no rescaling
i_p.addParamValue('no_gray', false, @islogical);
i_p.parse(imgpath, varargin{:});
imgpath = i_p.Results.imgpath;
range = i_p.Results.range;
no_gray = i_p.Results.no_gray;
I = imread(imgpath);
if (size(I, 3) > 1) && (~no_gray)
    I = rgb2gray(I);
end
if ~isa(I, 'double')
    I = double(I);
end
if ~isa(range, 'char')
    minval = min(min(I));
    maxval = max(max(I));
    low = range(1);     high = range(2);
    if ~no_gray
        % Grayscale (single channel) case
        I = (I - minval) / (maxval - minval);
        I = (I * (high - low)) + low;
    else
        % RGB (three channel) case
        w = size(I, 2); h = size(I, 1);
        I = (I - repmat(minval, h, w)) ./ (repmat(maxval - minval, h, w));
        I = (I * (high - low)) + low;
    end
end

end
