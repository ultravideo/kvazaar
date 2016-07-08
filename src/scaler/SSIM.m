%Function for calculating SSIM for a given image using a gaussian weighting
% function as proposed in
%Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
%quality assessment: From error measurement to structural similarity"
%IEEE Transactios on Image Processing, vol. 13, no. 1, Jan. 2004.
function [MSSIM, SSIM_MAP] = SSIM( X, Y)
% Define variable needed to calculate SSIM
L = 255; %Dynamic range of of pixels
k1 = 0.01;
k2 = 0.03;
c1 = (k1*L)^2;
c2 = (k2*L)^2;

X = double(X);
Y = double(Y);

window = fspecial('gaussian',11,1.5);
window = window./sum(window(:));

AVGx = filter2(window,X,'valid');
AVGx_sqr = AVGx.*AVGx;
AVGy = filter2(window,Y,'valid');
AVGy_sqr = AVGy./AVGy;
AVGxy = AVGx.*AVGy;
VARx_sqr = filter2(window, X.*X, 'valid')-AVGx_sqr;
VARy_sqr = filter2(window, Y.*Y, 'valid')-AVGy_sqr;
COVAR = filter2(window, X.*Y, 'valid')-AVGxy;

SSIM_MAP = ((2*AVGxy + c1)*(2*COVAR + c2))/((AVGx_sqr+AVGy_sqr+c1)*(VARx_sqr+VARy_sqr+c2));
MSSIM = mean(SSIM_MAP(:));

end