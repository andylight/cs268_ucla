function [M1, M2, F1, F2] = compute_matches(I1, P1, I2, P2, varargin)
%COMPUTE_MATCHES Estimates matching keypoints between I1 and I2.
%Input:
%   image I1
%   matrix P1: Nx2
%   image I2
%   matrix P2: Mx2
%Output:
%   [M1, M2, F1, F2] = compute_matches(...)
%       M1(i) = j if point p_i from I1 corresponds to point p_j from I2,
%       M1(i) = NaN if no correspondence is found
%       Fi(i) = NaN if no correspondence found, or an error metric
%% Arguments
i_p = inputParser;
i_p.addRequired('I1', @isnumeric);
i_p.addRequired('P1', @isnumeric);
i_p.addRequired('I2', @isnumeric);
i_p.addRequired('P2', @isnumeric);
i_p.addParamValue('T', 0.02, @isnumeric);        % Threshold for matching
i_p.addParamValue('w_win', 20, @isnumeric);     % Width of patch window
i_p.addParamValue('h_win', 20, @isnumeric);     % Height
i_p.addParamValue('method', 'ssd', @ischar);    % Patch similarity method
i_p.addParamValue('interactive', false, @islogical);
i_p.parse(I1, P1, I2, P2, varargin{:});
I1 = i_p.Results.I1;
P1 = i_p.Results.P1;
I2 = i_p.Results.I2;
P2 = i_p.Results.P2;
T = i_p.Results.T;
w_win = i_p.Results.w_win;  h_win = i_p.Results.h_win;
method = i_p.Results.method;
interactive = i_p.Results.interactive;
%%
w1 = size(I1, 2);   h1 = size(I1, 1);
w2 = size(I2, 2);   h2 = size(I2, 1);

M1 = nan([1, size(P1, 1)]);     F1 = nan([1, size(P1, 1)]);
M2 = nan([1, size(P2, 1)]);     F2 = nan([1, size(P2, 1)]);
if interactive == true
    figure;
end
for i=1:size(P1, 1)
    x_i = P1(i, 1);     y_i = P1(i, 2);
    P_i = get_img_patch(I1, x_i, y_i, w_win, h_win, 'fillval', nan);
    err_best = nan; j_best = nan;
    for j=1:size(P2, 1)
        x_j = P2(j, 1);     y_j = P2(j, 2);
        P_j = get_img_patch(I2, x_j, y_j, w_win, h_win, 'fillval', nan);
        if strcmp(method, 'ssd') == 1
            [err, P_diff] = imgdiff_L2(P_i, P_j);
        elseif strcmp(method, 'ncc')
            P_i_nonan = remove_border_nans(P_i);
            P_j_nonan = remove_border_nans(P_j);
            err = imgdiff_ncc(P_i_nonan, P_j_nonan);
        end
        %% START interactive
        if (interactive == true) && (err <= T)
            subplot(1, 3, 1);
            subimage(prepimage(P_i));
            hold on;
            plot(P1(1, :), P1(2, :), 'rx');
            hold off;
            title(sprintf('Patch i=%d', i));
            subplot(1, 3, 2);
            subimage(prepimage(P_j));
            hold on;
            plot(P2(1, :), P2(2, :), 'bx');
            hold off;
            title(sprintf('Patch j=%d', j));
            subplot(1, 3, 3);
            subimage(prepimage(P_diff));
            title(sprintf('Patch diff (err = %.4f)', err));
            if err <= T
                suptitle(sprintf('POTENTIAL MATCH (err=%.4f, T=%.4f)', err, T));
            else
                suptitle(sprintf('REJECT (err=%.4f, T=%.4f)', err, T));
            end
            pause;
        end
        %% END interactive
        if isnan(err_best) && err <= T
            err_best = err; j_best = j;
        elseif err < err_best
            err_best = err; j_best = j;
        end
    end
    if ~isnan(err_best)
        M1(i) = j_best;  M2(j_best) = i;
        F1(i) = err_best;    F2(j_best) = err_best;
    end
end
end

function Aout = remove_border_nans(A)
%REMOVE_BORDER_NANS Outputs a (potentially) smaller matrix AOUT that is
%the result of removing any rows/cols of A that are purely NaN.
Atmp0 = []; flag_foundval = false; cur_i = 1;
for i=1:size(A, 1)
    if ~all(isnan(A(i, :))) || (flag_foundval == true)
        flag_foundval = true;
        %Atmp0(cur_i, :) = A(i, :);
        Atmp0 = [Atmp0; A(i, :)];
        cur_i = cur_i + 1;
    end
end
Atmp1 = []; flag_foundval = false; cur_i = size(A, 1);
for i=size(Atmp0, 1):-1:1
    if ~all(isnan(Atmp0(i, :))) || (flag_foundval == true)
        flag_foundval = true;
        %Atmp1(cur_i, :) = A(i, :);
        %Atmp1 = [Atmp1; A(i, :)];
        Atmp1 = [Atmp0(i, :); Atmp1];
        cur_i = cur_i - 1;
    end
end
Atmp2 = []; flag_foundval = false; cur_j = 1;
for j=1:size(Atmp1, 2)
    if ~all(isnan(Atmp1(:, j))) || (flag_foundval == true)
        flag_foundval = true;
        %Atmp2(:, cur_j) = A(:, j);
        Atmp2 = [Atmp2, Atmp1(:, j)];
        cur_j = cur_j + 1;
    end
end
Atmp3 = []; flag_foundval = false; cur_j = size(A, 2);
for j=size(Atmp2, 2):-1:1
    if ~all(isnan(Atmp2(:, j))) || (flag_foundval == true)
        flag_foundval = true;
        %Atmp3(:, cur_j) = A(:, j);
        Atmp3 = [Atmp2(:, j), Atmp3];
        cur_j = cur_j + 1;
    end
end
Aout = Atmp3;
end