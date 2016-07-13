% Test SSIM of scaled image
% init mex
mex -g downScaler.c %-v COMPFLAGS='$COMPFLAGS /E' downScaler.c

%% Load an image and call downScaler

%rgb = imread('ngc6543a.jpg');
rgb = imread('peppers.png');
rgb = rgb(1:300,1:300,:);%250:300,275:325,:);
yuv = rgb2ycbcr(rgb);
factor = 2;
s = uint32(size(yuv(:,:,1))./factor);%uint32(size(yuv(275:300,300:325,1))./factor);

%% Downscale
[y,u,v] = downScaler( yuv(:,:,1), s, yuv(:,:,2), s, yuv(:,:,3), s);

%% Do a bad upsampling to get picture back to the same size

scaled = y;
scaled(:,:,2) = u;
scaled(:,:,3) = v;

scaled_rgb = ycbcr2rgb(scaled);

upscaled = zeros(size(rgb));
h = size(rgb,1);
w = size(rgb,2);
% Store SSIM for each layer
M = [];
I = zeros(3,1);
for n = 1:3
    upscaled(1:2:h, 1:2:h, n) = scaled(:,:,n);
    upscaled(2:2:h, 1:2:h, n) = scaled(:,:,n);
    upscaled(1:2:h, 2:2:h, n) = scaled(:,:,n);
    upscaled(2:2:h, 2:2:h, n) = scaled(:,:,n);
    
    [I(n,1),M(:,:,n)]=SSIM(yuv(:,:,n),upscaled(:,:,n));
end

upscaled_rgb = ycbcr2rgb(upscaled);

%% Display images
close all;
subplot(1,2,1);
imshow(rgb);
subplot(1,2,2);
imshow(upscaled);

figure;
subplot(1,3,1);
imshow(M(:,:,1));
subplot(1,3,2);
imshow(M(:,:,2));
subplot(1,3,3);
imshow(M(:,:,3));

disp(I);