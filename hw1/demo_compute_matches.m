function demo_compute_matches(varargin)
%DEMO_COMPUTE_MATCHES Summary of this function goes here
%   Detailed explanation goes here
i_p = inputParser;
i_p.addParamValue('imgsdir', 'imgs', @ischar);
i_p.addParamValue('ptsdir', 'ptsdata', @ischar);
i_p.addParamValue('ids', [1 2], @isnumeric);
i_p.addParamValue('T', 0.02, @isnumeric);
i_p.addParamValue('w_win', 20, @isnumeric);
i_p.addParamValue('h_win', 20, @isnumeric);
i_p.addParamValue('method', 'ssd', @ischar);
i_p.addParamValue('interactive', false, @islogical);
i_p.parse(varargin{:});
imgsdir = i_p.Results.imgsdir;
ptsdir = i_p.Results.ptsdir;
ids = i_p.Results.ids;
w_win = i_p.Results.w_win;      h_win = i_p.Results.h_win;
T = i_p.Results.T;
method = i_p.Results.method;
interactive = i_p.Results.interactive;

imgpaths = sort(getAllFiles(imgsdir));
imgpath1 = imgpaths{ids(1)};   imgpath2 = imgpaths{ids(2)};
disp(sprintf('Finding matches between images: %s %s', imgpath1, imgpath2));
datapath1 = get_data_path(imgpath1, ptsdir);
datapath2 = get_data_path(imgpath2, ptsdir);

I1 = imread_gray(imgpath1, 'range', [0 1]);
I2 = imread_gray(imgpath2, 'range', [0 1]);
imgs = {I1 I2};
P1 = load(datapath1);           P2 = load(datapath2);

[M1, M2, F1, F2] = compute_matches(I1, P1, I2, P2, ...
                           'w_win', w_win, 'h_win', h_win, ...
                           'T', T, ...
                           'method', method, ...
                           'interactive', interactive);

nb_matches = sum(~isnan(M1));
disp(sprintf('Found %d matches', nb_matches));
PtsMat1 = [];   PtsMat2 = [];
ErrsMat = [];
for i=1:length(M1)
    if ~isnan(M1(i))
        PtsMat1 = [PtsMat1; P1(i, :)];
        PtsMat2 = [PtsMat2; P2(M1(i), :)];
        ErrsMat = [ErrsMat F1(i)];
    end
end

if numel(PtsMat1) ~= 0 
    figure;
    subplot(2, 1, 1);
    showMatchedFeatures(I1, I2, PtsMat1, PtsMat2, 'montage');
    title('Putative point matches');
    legend('matchedPts1', 'matchedPts2');
    subplot(2, 1, 2);
    plot(1:length(ErrsMat), ErrsMat, 'bo');
    hold on;
    line([0 length(ErrsMat)+1], [T T], 'Color', 'red', 'LineStyle', '--');
    hold off;
    legend('Match Errors', sprintf('T=%.2f', T));
    axis([0 length(ErrsMat)+1 0.0 2*T]);
end

% G_21 := Transforms I2 -> I1
G_21 = recover_transform(PtsMat1, PtsMat2);
G_all = {eye(3), G_21};

% Lay out all I1, I2 on canvas and apply G to I2.
[imgdims_all, w_total, h_total] = get_all_imgdims(imgs);
[x0, y0] = get_origin(G_all);
Iall = zeros([h_total - y0 + 1, w_total - x0 + 1]);
for i=1:length(imgs)
    I = imgs{i};    G_i = G_all{i};
    [theta, tx, ty] = get_rigid_params(G_i);
    Trot = [[cos(theta) -sin(theta) 0];
            [sin(theta) cos(theta) 0];
            [0 0 1]];
    Irot = imwarp(I, affine2d(Trot'));
    w = size(Irot, 2); h = size(Irot, 1);
    i0 = ty - y0 + 1; i1 = i0 + h - 1;
    j0 = tx - x0 + 1; j1 = j0 + w - 1;
    Iall(i0:i1, j0:j1) = Irot;
end

figure;
imshow(Iall, []);

end

function [dims, w_total, h_total] = get_all_imgdims(imgs)
dims = {}; w_total = 0; h_total = 0;
for i=1:length(imgs)
    I = imgs{i};
    w = size(I, 2); h = size(I, 1);
    dims = {dims [w, h]};
    w_total = w_total + w;
    h_total = h_total + h;
end
end
function [theta, tx, ty] = get_rigid_params(G)
theta = acos(min(1, G(1,1)));   % Avoid taking acos(x) for x >= 1 (returns complex values)
tx = G(1, 3); ty = G(2, 3);
end
function [x0, y0] = get_origin(Gs)
x0 = 0; y0 = 0;
for i=1:length(Gs)
    G = Gs{i};
    [theta, tx, ty] = get_rigid_params(G);
    x0 = min(x0, tx);
    y0 = min(y0, ty);
end
end