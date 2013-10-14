function test_recover_transform()
%TEST_RECOVER_TRANSFORM Summary of this function goes here
%   Detailed explanation goes here
A = [[0 1];
     [1 0];
     [1 1];
     [-1 -2]
     [5 8];
     [11 12];
     [2 3]];
 
theta_opt = pi/4;
scale_opt = 2.0;
tx_opt = 10;
ty_opt = -4;

T_opt = [[scale_opt*cos(theta_opt) -scale_opt*sin(theta_opt) tx_opt];
         [scale_opt*sin(theta_opt) scale_opt*cos(theta_opt) ty_opt];
         [0 0 1]];
     
B = zeros([size(A, 1), 2]);
for i=1:size(A, 1)
    pt = A(i, :);
    pt_trans = T_opt * [pt 1]';
    B(i, :) = pt_trans(1:2);
end

% Add outlier to A, B
%A = [A; [4 4]];
%B = [B; [-2 -3]];

% Add gaussian noise to B
sigma = 0.1;
B = B + sigma*randn(size(B));

figure; subplot(1, 1, 1);
plot(A(:, 1), A(:, 2), 'ro'); hold on;
plot(B(:, 1), B(:, 2), 'bx');
legend('Orig. Points', 'Trans. Points');

% T transforms A to B
T = recover_transform(B, A);
theta = atan(T(2,1)/T(1,1));
scale = T(1,1) / cos(theta);
tx = T(1,3);
ty = T(2,3);

title(sprintf('Expect: theta=%.4f scale=%.4f tx=%.2f ty=%.2f\nGot: theta=%.4f scale=%.4f tx=%.2f ty=%.2f', ...
              theta_opt, scale_opt, tx_opt, ty_opt, theta, scale, tx, ty));

          
disp(sprintf('Expect: theta=%.4f scale=%.4f tx=%.2f ty=%.2f\nGot: theta=%.4f scale=%.4f tx=%.2f ty=%.2f', ...
              theta_opt, scale_opt, tx_opt, ty_opt, theta, scale, tx, ty));
          
B_out = zeros(size(B));
for i=1:size(A, 1)
    pt = [A(i, :) 1];
    pt_trans = T * pt';
    B_out(i, :) = pt_trans(1:2);
end

my_err = norm(sum(B - B_out), 2);
disp(sprintf('Total Error: %.4f', my_err));
end

