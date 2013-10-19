function G = stitch_graph(rootid, graph, matches, pts, history, varargin)
%STITCH_GRAPH Stitch together all images connected to image ROOT, and
%output the transformations into G.
%Input
%   int ROOTID
%   matrix GRAPH: (NxN) adjacency matrix
%   cell MATCHES: {NxN} cell where
%       MATCHES{i,j} = {Mi, Mj, Fi, Fj}
%   cell PTS: {1xN} cell where
%       PTS{i} = (Mx2) matrix
%   containers.Map HISTORY: (int32 -> logical)
%       Stores node indices that have already been visited.
%Output
%   cell G: (1xN) cell, where
%       {{int imgid, matrix T}, ...};
history(int32(rootid)) = true;
G = {{rootid, eye(3)}};
neighbors = get_neighbors(graph, rootid);
if numel(neighbors) == 0
    return;
end
for i=1:numel(neighbors)
    neighborid = neighbors{i};
    if history.isKey(neighborid)
        continue;
    end
    G_sub = stitch_graph(neighborid, graph, matches, pts, history);
    matcell = matches{rootid, neighborid};
    M1 = matcell{1}; M2 = matcell{2};
    pts1 = pts{rootid}; pts2 = pts{neighborid};
    [matpts1, matpts2] = get_matched_pts(pts1, pts2, M1, M2);
    T = recover_transform(matpts1, matpts2);
    err = compute_err(matpts1, matpts2, T);
    %disp(sprintf('err is: %.6f', err));
    for ii=1:numel(G_sub)
        imgid = G_sub{ii}{1};   T_sub = G_sub{ii}{2};
        G{numel(G)+1} = {imgid, T * T_sub};
    end
end
end

function [matpts1, matpts2] = get_matched_pts(pts1, pts2, M1, M2)
%GET_MATCHED_PTS Given point correspondences M1,M2, output the matching
%points.
nb_matches = sum(~isnan(M1));
matpts1 = zeros([nb_matches, 2]);   matpts2 = zeros([nb_matches, 2]);
cnt = 1;
for i=1:length(M1)
    val = M1(i);
    if ~isnan(val)
        pt1 = pts1(i, :);
        pt2 = pts2(val, :);
        matpts1(cnt, :) = pt1;
        matpts2(cnt, :) = pt2;
        cnt = cnt + 1;
    end
end
end

function err = compute_err(pts1, pts2, T)
err = 0;
for i=1:size(pts1, 1)
    pt1 = pts1(i, :); pt2 = pts2(i, :);
    pt1 = [pt1 1];  pt2 = [pt2 1];
    err = err + sum(((T * pt2') - pt1').^2);
end
end