%% Salt and pepper noise SSIM test
%
rgb = imread('ngc6543a.jpg');
ref = rgb;
rate = 0.001;
s = size(rgb(:,:,1));

% Add 1 and 0 values based at random

inds = repmat(rand(s) < rate,1,1,3);
rgb(inds) = 255;
inds = repmat(rand(s) < rate,1,1,3);
rgb(inds) = 0;

%% Show results
imshow(rgb);