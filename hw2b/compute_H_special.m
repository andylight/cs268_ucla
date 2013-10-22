function H = compute_H_special(p,q)
%COMPUTE_H_SPECIAL Summary of this function goes here
%   p,q := 3xN matrix of 2d homogenous points
n = size(p, 2);
chi = zeros([9, 3*n]);
for j=1:n
    pj = p(:, j);
    qj = q(:, j);
    qj_hat = make_hat(qj);
    aj = kron(qj_hat, pj);
    i_start = 3*(j-1) + 1;
    i_end = i_start + 2;
    chi(:, i_start:i_end) = aj;
end
chi = chi';

X1 = p';
X2 = q';
chi2 = make_chi(X1, X2);
chi = chi2;

[U, S, Vt] = svd(chi);
V = Vt';
HL = V(9,:);
H = reshape(HL, [3 3]); % rearranges HL into columns of H
end

function chi = make_chi(X1, X2)
% X1,X2 are Nx3 matrix of 2d points in homogeneous coordinates
chi = [];
n = size(X1, 1);
for i=1:n
    x1 = X1(i,:)';
    x2 = X2(i,:)';
    x2_hat = make_hat(x2);
    a = kron(x1, x2_hat);
    chi = [chi; a'];
end
end

function xhat = make_hat(x)
xhat = zeros([3,3]);
xhat(1,2) = -x(3);
xhat(1,3) = x(2);
xhat(2,1) = x(3);
xhat(2,3) = -x(1);
xhat(3,1) = -x(2);
xhat(3,2) = x(1);
end
