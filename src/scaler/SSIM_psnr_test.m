%% Test the effect of quality loss from sequential down- and upscaling to PSNR/MSSIM
% MEX compilation
mex -g downScaler.c

%% Load test image

rgb = imread('peppers.png');%'ngc6543a.jpg');
yuv = rgb2ycbcr(rgb);

%% Calculate metrics
limit = 0.05; %How low we try to scale at most ( size.*limit )
num = 100; %Number of test points

s = size(rgb(:,:,1));
A = [linspace(s(1)*limit,s(1),num) ; linspace(s(2)*limit,s(2),num)];
ratio = A(1,:)./s(1); %; A(2,:)./s(2)];
s = uint32(s);
A = uint32(A);

SI = zeros(3,size(A,2));
FSI = SI;
PI = SI;
ind = 1;
for r = A
    % Downscale then upscale
    [y, u, v] = downScaler( yuv(:,:,1), r', yuv(:,:,2), r', yuv(:,:,3), r' );
    [y, u, v] = downScaler( y, s, u, s, v, s );
    
    % Calculate PSNR/SSIM
    PI(1,ind) = ssim(y,yuv(:,:,1));
    [SI(1,ind), ~] = SSIM(y, yuv(:,:,1));
    [FSI(1,ind), ~] = FastSSIM(y, yuv(:,:,1));
    PI(2,ind) = ssim(u,yuv(:,:,2));
    [SI(2,ind), ~] = SSIM(u, yuv(:,:,2));
    [FSI(2,ind), ~] = FastSSIM(u, yuv(:,:,2));
    PI(3,ind) = ssim(v,yuv(:,:,3));
    [SI(3,ind), ~] = SSIM(v, yuv(:,:,3));
    [FSI(3,ind), ~] = FastSSIM(v, yuv(:,:,3));
    ind = ind+1;
end

%% Plot results
close all;

plot(ratio,PI(1,:));
hold all;
plot(ratio,PI(2,:));
plot(ratio,PI(3,:));
title('PSNR');
legend('Y','U','V');

figure;

plot(ratio,SI(1,:));
hold all;
plot(ratio,SI(2,:));
plot(ratio,SI(3,:));
title('MSSIM');
legend('Y','U','V');

figure;

plot(ratio,FSI(1,:));
hold all;
plot(ratio,FSI(2,:));
plot(ratio,FSI(3,:));
title('Fast-MSSIM');
legend('Y','U','V');