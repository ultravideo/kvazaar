%% Run test on downScaler (4:2:0 chroma format + others)
% init mex
mex -g downScaler.c %-v COMPFLAGS='$COMPFLAGS /E' downScaler.c

%% Load an image and call downScaler

rgb = imread('ngc6543a.jpg');

%rgb = rgb(1:602,1:601,:);%250:300,275:325,:);
yuv = rgb2ycbcr(rgb);
factor = 2;
s = uint32([403 502]);%size(yuv(:,:,1))./factor);%uint32(size(yuv(275:300,300:325,1))./factor);
chroma_inds_h = round(linspace(1,size(yuv,1),size(yuv,1)/2));
chroma_inds_w = round(linspace(1,size(yuv,2),size(yuv,2)/2));
%% Downscale
[y,u,v] = downScaler( yuv(:,:,1), s, yuv(chroma_inds_h,chroma_inds_w,2), idivide(s,2,'floor'), yuv(chroma_inds_h,chroma_inds_w,3), idivide(s,2,'floor'));

%% Display images
subplot(1,3,1);
imshow(y);
subplot(1,3,2);
imshow(u);
subplot(1,3,3);
imshow(v);

%% Test 4:2:2
d = uint32([1 2]);
[y,u,v] = downScaler( yuv(:,:,1), s, yuv(:,chroma_inds_w,2), idivide(s,d,'floor'), yuv(:,chroma_inds_w,3), idivide(s,d,'floor'));

%% Display images
figure;
subplot(1,3,1);
imshow(y);
subplot(1,3,2);
imshow(u);
subplot(1,3,3);
imshow(v);

%% Test 4:0:0
[y,u,v] = downScaler( yuv(:,:,1), s, uint8([]), [0 0], uint8([]), [0 0]);

%%
figure;
imshow(y);