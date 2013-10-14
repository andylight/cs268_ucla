function wtf()
%WTF Summary of this function goes here
%   Detailed explanation goes here

I1 = imread_gray('imgs/I_0.png', 'range', [0 1], 'no_gray', true);
I2 = imread_gray('imgs/I_1.png', 'range', [0 1], 'no_gray', true);
I3 = imread_gray('imgs/I_2.png', 'range', [0 1], 'no_gray', true);

P1 = load('pts/I_0.pts');
P2 = load('pts/I_1.pts');
P3 = load('pts/I_2.pts');

G_11 = eye(3);
G_12 =  [[0.5284   -0.6986  195.8645];
         [0.6986    0.5284 -175.8208];
         [0         0    1.0000]];
G_12_notrans = [[G_12(1,1), G_12(1,2), 0];
                [G_12(2,1), G_12(2,2), 0]; [0 0 1]];
% theta: 52.89 deg scale=0.8759
G_13 =  [[-0.2822 -0.7045 779.8810];
         [0.7045 -0.2822 360.7269];
         [0 0 1]];
G_13_notrans = [[G_13(1,1), G_13(1,2), 0]; [G_13(2,1), G_13(2,2), 0]; [0 0 1]];
% theta: -68.16 deg scale=-0.759
M_A_11 = [1 2 3 nan 5];     M_A_11 = M_A_11(~isnan(M_A_11));
M_A_12 = [1 2 3 nan 5];     M_A_12 = M_A_12(~isnan(M_A_12));

M_B_11 = [1 2 3 4 nan];     M_B_11 = M_B_11(~isnan(M_B_11));
M_B_13 = [1 2 3 4];         M_B_13 = M_B_13(~isnan(M_B_13));
     
PtsMat_A_11 = P1(M_A_11, :);
PtsMat_A_12 = P2(M_A_12, :);

PtsMat_B_11 = P1(M_B_11, :);
PtsMat_B_13 = P3(M_B_13, :);

IA_1 = I1;
IA_2 = imwarp(I2, affine2d(G_12'));
IB_1 = I1;
IB_3 = imwarp(I3, affine2d(G_13'));

PtsMat_A_12_trans = zeros(size(PtsMat_A_12));
for i=1:size(PtsMat_A_12, 1)
    pt = PtsMat_A_12(i, :);
    pt_t = G_12 * [pt 1]';
    pt_t(1) = pt_t(1) + 297 + 195;
    pt_t(2) = pt_t(2) + 81 - 175;
    PtsMat_A_12_trans(i, :) = pt_t(1:2);
end

PtsMat_B_13_trans = zeros(size(PtsMat_B_13));
for i=1:size(PtsMat_B_13, 1)
    pt = PtsMat_B_13(i, :);
    pt_t = G_13 * [pt 1]';
    PtsMat_B_13_trans(i, :) = pt_t(1:2);
end

PtsMat_A_11 = [PtsMat_A_11; 1 1];
PtsMat_A_12_trans = [PtsMat_A_12_trans; 1 1];

I2_ctr = [size(I2, 2) / 2, size(I2, 1) / 2];
I3_ctr = [size(I3, 2) / 2, size(I3, 1) / 2];

I2_ctr_trans = G_12 * [I2_ctr 1]';  I2_ctr_trans(1:2);
I3_ctr_trans = G_13 * [I3_ctr 1]';  I3_ctr_trans(1:2);

I2_xdiff = I2_ctr_trans(1) - I2_ctr(1);
I2_ydiff = I2_ctr_trans(2) - I2_ctr(2);
I2_realtrans = [G_12(1,3) + I2_xdiff, G_12(2,3) + I2_ydiff];

I3_xdiff = I3_ctr_trans(1) - I3_ctr(1);
I3_ydiff = I3_ctr_trans(2) - I3_ctr(2);
I3_realtrans = [G_13(1,3) + I3_xdiff, G_13(2,3) + I3_ydiff];

disp(sprintf('I2 origin mvmt: tx=%.4f ty=%.4f', I2_xdiff, I2_ydiff));
I2_realtrans
disp(sprintf('I3 origin mvmt: tx=%.4f ty=%.4f', I3_xdiff, I3_ydiff));
I3_realtrans

figure;
showMatchedFeatures(IA_1, IA_2, PtsMat_A_11, PtsMat_A_12_trans, 'montage');
title(sprintf('I1, I2. tx=%.4f ty=%.4f', G_12(1,3), G_12(2,3)));
legend('I1_pts', 'I2_pts');
figure;
showMatchedFeatures(IB_1, IB_3, PtsMat_B_11, PtsMat_B_13_trans, 'montage');
title(sprintf('I1, I3. tx=%.4f ty=%.4f', G_13(1,3), G_13(2,3)));
legend('I1_pts', 'I3_pts');

RA = imref2d(size(I2), [0 1000], [0 1000]);

[B, RB] = imwarp(I2, affine2d(G_12'));
figure;
imshow(B, RB);

end

