function [M1, M2, varargout] = compute_corr_simple(P, Q, varargin)
%COMPUTE_CORR_SIMPLE Computes correspondences between points on image
%I1 and points on image I2.
%Input:
%   matrix P: Nx2
%   matrix Q: Mx2
%Output:
%   [M1, M2] = compute_correspondences(P, Q)
%       M1(i) = j if point P(i) corresponds to point Q(j)
%             = NaN if P(i) does not correspond with Q(j)
%       M2(i) = j if point Q(i) corresponds to point P(j)
%             = NaN if Q(i) does not correspond with P(j)
%   [M1, M2, T] = compute_correspondences(P, Q)
%       T is the computed transformation that maps Q to P.
N = 5;    % Min. number of data required to fit model
K = 10000; % Number of iterations
best_model = 0;
best_error = inf;
best_idxs1 = 0;
best_idxs2 = 0;

for cur_iter=1:K
    % maybe_inliers := {P*, Q*} where P*,Q* are Nx2 points
    %   idxs1 := 1xN vector of indexes into P
    [maybe_inliers, idxs1, idxs2] = randselect(P, Q, N);
    P_inliers = maybe_inliers{1};
    Q_inliers = maybe_inliers{2};
    % T aligns points Q to P
    [T, model_err] = estimate_model(P_inliers, Q_inliers);
    % If this is the best model, keep it
    if model_err <= best_error
        best_model = T;
        best_idxs1 = idxs1;
        best_idxs2 = idxs2;
        best_error = model_err;
    end
end
M1 = nan([1, size(P, 1)]);
M2 = nan([1, size(Q, 1)]);
for i=1:size(best_idxs2, 2)
    M1(best_idxs1(i)) = best_idxs2(i);
    M2(best_idxs2(i)) = best_idxs1(i);
end
if nargout >= 3
    varargout{1} = best_model;
end
if nargout >= 4
    varargout{2} = best_error;
end
end

function [pts, idxs1, idxs2] = randselect(P, Q, N)
%RANDSELECT Randomly select N points between P and Q
%   [pts idxs1 idxs2] = randselect(P, Q, N)
%       PTS := {P*, Q*} where each P* is an Nx2 matrix of points
%       IDXS1 := 1xN vector containing the (row) indices from P used in P*.
N_min = min(N, size(P, 1));     % Don't over-choose points
idxs1 = randperm(size(P, 1));
idxs2 = randperm(size(Q, 1));
idxs1 = idxs1(1:N_min);
idxs2 = idxs2(1:N_min);
%idxs1 = sort(idxs1(1:N_min));
%idxs2 = sort(idxs2(1:N_min));
pts = {P(idxs1, :), Q(idxs2, :)};
end

function [G, err] = estimate_model(P1, P2)
%RECOVER_TRANSFORM Computes transformation that aligns points in P2 to P1.
% P1, P2 are Nx2 matrices.
if size(P1, 1) <= 1
    % Don't have enough points to estimate rigid transform. Do simple
    % translation estimate instead as a fallback.
    A = zeros([2*size(P2, 1), 3]);
    for i=1:size(P2, 1)
        pt = P2(i, :);
        x = pt(1); y = pt(2);
        idx1 = (i*2) - 1;
        idx2 = idx1 + 1;
        A(idx1, :) = [1 0 x];
        A(idx2, :) = [0 1 y];
    end
    b = reshape(P1', [size(P1, 1) * size(P1, 2), 1]);
    % Constrain that the last variable in p is 1, i.e.: p(3) == 1
    A = [A; 0 0 1];  b = [b; 1];
    p = A\b;
    tx = p(1); ty = p(2);
    G = [[1 0 tx];
         [0 1 ty];
         [0 0 1]];
    p_eval = [tx, ty, 1];
    err = norm(A*p_eval' - b, 2);
else
    A = zeros([2*size(P2, 1), 4]);
    for i=1:size(P2, 1)
        pt = P2(i, :);
        x = pt(1);  y = pt(2);
        idx1 = (i*2) - 1;
        idx2 = idx1 + 1;
        A(idx1, :) = [x -y 1 0];
        A(idx2, :) = [y x 0 1];
    end

    b = reshape(P1', [size(P1, 1) * size(P1, 2), 1]);

    p = A\b;
    alpha = p(1); beta = p(2);
    tx = p(3);  ty = p(4);
    G = [[alpha -beta tx];
         [beta alpha ty];
         [0 0 1]];
    err = norm(A*p - b, 2);
end
end

function out = transform_points(T, Q)
%TRANSFORM_POINTS Transforms points Q with transformation T
%   out = transform_points(T, Q)
%       T is a 3x3 affine matrix. Q is an Nx2 matrix of points.
%       OUT is a Nx2 matrix of transformed points.
Q_homo = pts2homo(Q);   % Nx3 matrix of points
Q_trans = T * Q_homo';  % 3xN matrix of transformed points
out = Q_trans(1:2, :);
out = out';
end

function out = pts2homo(Q)
%PTS2HOMO Converts points into homogenous coords.
out = zeros([size(Q, 1), size(Q, 2) + 1]);
for i=1:size(Q, 1)
    pt = Q(i, :);
    out(i,:) = [pt 1];
end
end

function err = compute_error(T, P, Q)
%COMPUTE_ERROR Computes the error of aligning points Q to P with transformation
%T.
err = 0;
Q_align = transform_points(T, Q);
for i=1:size(Q_align, 1)
    q_i = Q_align(i, :);
    [pclose_i, dist] = find_closest(P, q_i);
    err = err + dist;
end
end

function [p_idx, dist] = find_closest(P, q)
%FIND_CLOSEST Finds the point in P that is closest to point Q.
p_idx = nan;
dist = inf;
for i=1:size(P, 1)
    p_i = P(i, :);
    curdist = norm(p_i - q, 2);
    if curdist <= dist
        p_idx = i;
        dist = curdist;
    end
end
end
