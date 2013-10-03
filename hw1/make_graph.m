function [G, matches] = make_graph(imgs, pts, varargin)
%MAKE_GRAPH Constructs graph that encodes the spatial relationships between
%all images. Nodes are images. An edge is present between I_i, I_j if at
%least one point in I_i corresponds to a point in I_j.
%Input
%   cell IMGS: (1xN) cell array
%   cell PTS: (1xN) cell array, each containing (Mx2) matrices
% MATCHES: (NxN) cell array
%   MATCHES{i,j} = {Mi, Mj, Fi, Fj}
i_p = inputParser; 
i_p.addRequired('imgs');
i_p.addRequired('pts');
i_p.addParamValue('T', 0.05, @isnumeric);
i_p.addParamValue('method', 'ssd', @ischar);
i_p.parse(imgs, pts, varargin{:});
imgs = i_p.Results.imgs;
pts = i_p.Results.pts;
T = i_p.Results.T;
method = i_p.Results.method;
matches = cell([numel(imgs), numel(imgs)]);
G = zeros([numel(imgs), numel(imgs)]);
for i=1:length(imgs)
    I_i = imgs{i};
    pts_i = pts{i};
    for j=i+1:length(imgs)
        I_j = imgs{j};
        pts_j = pts{j};
        [M1, M2, F1, F2] = compute_matches(I_i, pts_i, I_j, pts_j, 'T', T, 'method', method);
        nb_matches = sum(~isnan(M1));
        matches{i,j} = {M1, M2, F1, F2};
        matches{j,i} = {M2, M1, F2, F1};
        if nb_matches >= 1
            G(i, j) = 1;
            G(j, i) = 1;
        end
    end
end
end
