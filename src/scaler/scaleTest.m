%% Run test on downScaler
% init mex
mex -g downScaler.c %-v COMPFLAGS='$COMPFLAGS /E' downScaler.c

%% Load an image and call downScaler

rgb = imread('saturn.png');%imread('ngc6543a.jpg');

%rgb = rgb(1:602,1:601,:);%250:300,275:325,:);
yuv = rgb2ycbcr(rgb);
factor = 2;
s = uint32([1000 1000]);%size(yuv(:,:,1))./factor);%uint32(size(yuv(275:300,300:325,1))./factor);

%% Downscale
[y,u,v] = downScaler( yuv(:,:,1), s, yuv(:,:,2), s, yuv(:,:,3), s);

%% Display images

scaled = y;
scaled(:,:,2) = u;
scaled(:,:,3) = v;

scaled_rgb = ycbcr2rgb(scaled);

subplot(1,2,1);
imshow(rgb);
subplot(1,2,2);
imshow(scaled_rgb);