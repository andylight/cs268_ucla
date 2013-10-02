function test_imgpaste()
%TEST_IMGPASTE Summary of this function goes here
%   Detailed explanation goes here
I = imread_gray('../../corgi_dapper.png', 'range', [0 1], 'no_gray', true);

x1 = 100; y1 = 125;
x2 = 140; y2 = 150;
AMT = 0.3;

patch1 = I(y1:y2, x1:x2, :);
patch1 = patch1 + AMT;
patch1(patch1 > 1) = 1; patch1(patch1 < 0) = 0;

Ipaste_overwrite = imgpaste(I, patch1, x1, y1, 'method', 'overwrite');
Ipaste_average = imgpaste(I, patch1, x1, y1, 'method', 'average');
Ipaste_meanshift = imgpaste(I, patch1, x1, y1, 'method', 'meanshift');

figure;
subplot(2, 3, 1);
imshow(I, []);
title('I');
subplot(2, 3, 2);
imshow(patch1, []);
title(sprintf('patch adjusted by: %d', AMT));
subplot(2, 3, 4);
imshow(prepimage(Ipaste_overwrite), []);
title('overwrite');
subplot(2, 3, 5);
imshow(prepimage(Ipaste_average), []);
title('average');
subplot(2, 3, 6);
imshow(prepimage(Ipaste_meanshift), []);
title('Mean Shift');

% Grayscale test
Igray = imread_gray('../../corgi_dapper.png', 'range', [0 1]);
patch1 = Igray(y1:y2, x1:x2);
patch1 = patch1 + AMT;
patch1(patch1 > 1) = 1; patch1(patch1 < 0) = 0;

Ipaste_overwrite = imgpaste(Igray, patch1, x1, y1, 'method', 'overwrite');
Ipaste_average = imgpaste(Igray, patch1, x1, y1, 'method', 'average');
Ipaste_meanshift = imgpaste(Igray, patch1, x1, y1, 'method', 'meanshift');

figure;
subplot(2, 3, 1);
imshow(Igray, []);
title('Igray');
subplot(2, 3, 2);
imshow(patch1, []);
title(sprintf('patch adjusted by: %d', AMT));
subplot(2, 3, 4);
imshow(prepimage(Ipaste_overwrite), []);
title('overwrite');
subplot(2, 3, 5);
imshow(prepimage(Ipaste_average), []);
title('average');
subplot(2, 3, 6);
imshow(prepimage(Ipaste_meanshift), []);
title('Mean Shift');

end

