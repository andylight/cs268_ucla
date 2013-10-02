function G = recover_transform(P1, P2)
%RECOVER_TRANSFORM Computes transformation that aligns points in P2 to P1.
% P1, P2 are Nx2 matrices.
if size(P1, 1) < 4
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
    return;
end
A = zeros([2*size(P2, 1), 4]);
for i=1:size(P2, 1)
    pt = P2(i, :);
    x = pt(1);  y = pt(2);
    idx1 = (i*2) - 1;
    idx2 = idx1 + 1;
    A(idx1, :) = [-y x 1 0];
    A(idx2, :) = [x y 0 1];
end

b = reshape(P1', [size(P1, 1) * size(P1, 2), 1]);

p = A\b;
s = p(1);   c = p(2);
tx = p(3);  ty = p(4);
% Enforce that c <= 1 to avoid bogus acos(theta) outputs.
G = [[min(1, c) -s tx];
     [s min(1, c) ty];
     [0 0 1]];
end

