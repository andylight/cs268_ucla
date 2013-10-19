function [theta, sc, tx, ty] = get_trans_params(G)
%GET_TRANS_PARAMS Returns the parameters of this transformation matrix,
%assuming SIMILARITY.
theta = atan(G(2,1)/G(1,1));
sc = G(1,1) / cos(theta);

G_rot_scale = [[G(1,1), G(1,2), 0];
               [G(2,1), G(2,2), 0];
               [0 0 1]];
G_trans = inv(G_rot_scale) * G;
tx = G_trans(1,3);
ty = G_trans(2,3);
end
