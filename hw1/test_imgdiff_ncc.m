function test_imgdiff_ncc()
%TEST_IMGDIFF_NCC Summary of this function goes here
%   Detailed explanation goes here
imgpath1 = 'imgs/I_02.png';
imgpath2 = 'I_02_sub_x20_y25.png';

I1 = imread_gray(imgpath1, 'range', [0, 1]);
I2 = imread_gray(imgpath2, 'range', [0, 1]);

[err, F] = imgdiff_ncc(I1, I2, 'interactive', true);

[err2, F2] = imgdiff_ncc(imread_gray('corgi1.png', 'range', [0, 1]), ...
                         imread_gray('corgi2.png', 'range', [0, 1]), ...
                         'interactive', true);
[err2, F2] = imgdiff_ncc(imread_gray('corgi1.png', 'range', [0, 1]), ...
                         imread_gray('corgi4.png', 'range', [0, 1]), ...
                         'interactive', true);
                 

end

