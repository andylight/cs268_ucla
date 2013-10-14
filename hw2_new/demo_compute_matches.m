function demo_compute_matches(varargin)
%DEMO_COMPUTE_MATCHES Summary of this function goes here
%   Detailed explanation goes here
i_p = inputParser;
i_p.addParamValue('imgsdir', 'imgs', @ischar);
i_p.addParamValue('ptsdir', 'pts', @ischar);
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

%[M1, M2, F1, F2] = compute_matches(I1, P1, I2, P2, ...
%                           'w_win', w_win, 'h_win', h_win, ...
%                           'T', T, ...
%                           'method', method, ...
%                           'interactive', interactive);
[M1, M2, Tmat, err] = compute_corr_simple(P1, P2);
F1 = 0; F2 = 0;
[th, sc, tx, ty] = get_trans_params(Tmat);
disp(sprintf('Got: theta=%.4f scale=%.4f tx=%.2f ty=%.2f', rad2deg(th), sc, tx, ty));
disp(Tmat);
Trotscale = [[Tmat(1,1), Tmat(1,2), 0];
             [Tmat(2,1), Tmat(2,2), 0];
             [0 0 1]];
disp(Trotscale*[tx ty 1]');

nb_matches = sum(~isnan(M1));
disp(sprintf('Found %d matches', nb_matches));
PtsMat1 = [];   PtsMat2 = [];
ErrsMat = [];
for i=1:length(M1)
    if ~isnan(M1(i))
        PtsMat1 = [PtsMat1; P1(i, :)];
        PtsMat2 = [PtsMat2; P2(M1(i), :)];
        %ErrsMat = [ErrsMat F1(i)];
    end
end
M1
M2
if numel(PtsMat1) ~= 0 
    figure;
    subplot(2, 2, 1);
    showMatchedFeatures(I1, I2, PtsMat1, PtsMat2, 'montage');
    title('Putative point matches');
    legend('matchedPts1', 'matchedPts2');
    subplot(2, 2, 2);
    plot(1:length(ErrsMat), ErrsMat, 'bo');
    hold on;
    line([0 length(ErrsMat)+1], [T T], 'Color', 'red', 'LineStyle', '--');
    hold off;
    legend('Match Errors', sprintf('T=%.2f', T));
    axis([0 length(ErrsMat)+1 0.0 2*T]);
    subplot(2,2,3);
    PtsMat1_computed = zeros(size(PtsMat2));
    for i=1:size(PtsMat2, 1)
        pt2 = PtsMat2(i, :);
        pt1_comp = Tmat * [pt2 1]';
        PtsMat1_computed(i, :) = pt1_comp(1:2);
    end
    showMatchedFeatures(I1, I2, PtsMat1_computed, PtsMat2, 'montage');
    title('Computed Transformation');
    legend('pts1_computed', 'pts2');
    err = norm(PtsMat1_computed - PtsMat1, 2);
    disp(sprintf('l2 norm: %f', err));
end

% G_21 := Transforms I2 -> I1
G_21 = recover_transform(PtsMat1, PtsMat2);
G_all = {eye(3), G_21};

% Lay out all I1, I2 on canvas and apply G to I2.
[imgdims_all, w_total, h_total] = get_all_imgdims(imgs);
w_total = 1000; h_total = 1000;
[x0, y0] = get_origin(G_all);
x0 = -200;
y0 = -400;
Iall = zeros([h_total - y0 + 1, w_total - x0 + 1]);
for i=1:length(imgs)
    I = imgs{i};    G_i = G_all{i};
    [theta, sc, tx, ty] = get_trans_params(G_i);
    %Trot = [[cos(theta) -sin(theta) 0];
    %        [sin(theta) cos(theta) 0];
    %        [0 0 1]];
    Trot = [[G_i(1,1), G_i(1,2), 0];
            [G_i(2,1), G_i(2,2), 0];
            [0 0 1]];
    Irot = imwarp(I, affine2d(Trot'));
    w = size(Irot, 2); h = size(Irot, 1);
    %i0 = ty - y0 + 1; i1 = i0 + h - 1;
    %j0 = tx - x0 + 1; j1 = j0 + w - 1;
    new_xy = [G_i(1,3), G_i(2,3)];
    i0 = new_xy(2)-y0+1; i1 = i0 + h - 1;
    j0 = new_xy(1)-x0+1; j1 = j0 + w - 1;
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
function [theta, sc, tx, ty] = get_trans_params(G)
%GET_TRANS_PARAMS Returns the parameters of this transformation matrix.
theta = atan(G(2,1)/G(1,1));
sc = G(1,1) / cos(theta);

G_rot_scale = [[G(1,1), G(1,2), 0];
               [G(2,1), G(2,2), 0];
               [0 0 1]];
G_trans = inv(G_rot_scale) * G;
tx = G_trans(1,3);
ty = G_trans(2,3);
end
function [x0, y0] = get_origin(Gs)
x0 = 0; y0 = 0;
for i=1:length(Gs)
    G_mat = Gs{i};
    [theta, sc, tx, ty] = get_trans_params(G_mat);
    G_rot_scale = [[G_mat(1,1), G_mat(1,2), 0];
                   [G_mat(2,1), G_mat(2,2), 0];
                   [0 0 1]];
    offset = G_rot_scale*[1 1 1]';
    x0_off = offset(1) - 1;
    y0_off = offset(2) - 1;
    x0 = min(x0, abs(sc)*tx + x0_off);
    y0 = min(y0, abs(sc)*ty + y0_off);
end
end