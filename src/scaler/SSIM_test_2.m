% Test SSIM of downscaled then upscaled image
% init mex
mex -g downScaler.c %-v COMPFLAGS='$COMPFLAGS /E' downScaler.c

%% Load an image and call downScaler

rgb = imread('ngc6543a.jpg');

rgb = rgb(1:600,1:600,:);%250:300,275:325,:);
yuv = rgb2ycbcr(rgb);
factor = 2;
s1 = uint32([100 100]);
s2 = uint32([600 600]);%size(yuv(:,:,1))./factor);%uint32(size(yuv(275:300,300:325,1))./factor);

%% Downscale
[y1,u1,v1] = downScaler( yuv(:,:,1), s1, yuv(:,:,2), s1, yuv(:,:,3), s1);
[y2,u2,v2] = downScaler( y1, s2, u1, s2, v1, s2 );

%% Calculate SSIM

scaled = y1;
scaled(:,:,2) = u1;
scaled(:,:,3) = v1;

scaled2 = y2;
scaled2(:,:,2) = u2;
scaled2(:,:,3) = v2;

scaled_rgb = ycbcr2rgb(scaled);
scaled_rgb2 = ycbcr2rgb(scaled2);

% Store SSIM for each layer
M = [];
I = zeros(3,1);
PSNR = zeros(3,1);
for n = 1:3
    PSNR(n) = psnr(scaled2(:,:,n),yuv(:,:,n));
    [I(n,1),M(:,:,n)]=SSIM(yuv(:,:,n),scaled2(:,:,n));
end

%% Display images
close all;
subplot(1,2,1);
imshow(rgb);
subplot(1,2,2);
imshow(scaled_rgb2);

figure;
subplot(1,4,1);
imshow(M(:,:,1));
title(['MSSIM= ' num2str(I(1))]);
subplot(1,4,2);
imshow(M(:,:,2));
title(['MSSIM= ' num2str(I(2))]);
subplot(1,4,3);
imshow(M(:,:,3));
title(['MSSIM= ' num2str(I(3))]);
subplot(1,4,4);
imshow(M);

disp('MSSIM');
disp(I);
disp('Psnr');
disp(PSNR);