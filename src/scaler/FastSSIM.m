%Function for calculating Fast SSIM for a given image using a gaussian weighting
% function as proposed in
%Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
%quality assessment: From error measurement to structural similarity"
%IEEE Transactios on Image Processing, vol. 13, no. 1, Jan. 2004.
% and
%Chen, M. & Bovik, A.C. J Real-Time Image Proc (2011) 6: 281. doi:10.1007/s11554-010-0170-9
function [MSSIM, SSIM_MAP] = FastSSIM( X, Y)
% Define variable needed to calculate SSIM
L = 255; %Dynamic range of of pixels
k1 = 0.01;
k2 = 0.03;
c1 = (k1*L)^2;
c2 = (k2*L)^2;

X = double(X);
Y = double(Y);

%Integer aproximation for a 8x8 gaussian window
window = [0 0 0 1 1 0 0 0
          0 0 1 2 2 1 0 0 
          0 1 2 4 4 2 1 0
          1 2 4 8 8 4 2 1 
          1 2 4 8 8 4 2 1 
          0 1 2 4 4 2 1 0
          0 0 1 2 2 1 0 0 
          0 0 0 1 1 0 0 0];
%fspecial('gaussian',11,1.5);
%window = window/sum(window(:));

AVGx = integralMean(X);%filter2(window,X,'valid');
AVGx_sqr = AVGx.*AVGx;
AVGy = integralMean(Y);%filter2(window,Y,'valid');
AVGy_sqr = AVGy.*AVGy;
AVGxy = AVGx.*AVGy;
%VARx_sqr = filter2(window, X.*X, 'valid')-AVGx_sqr;
%VARy_sqr = filter2(window, Y.*Y, 'valid')-AVGy_sqr;
%COVAR = filter2(window, X.*Y, 'valid')-AVGxy;
GX = abs(gradMag(X));
GY = abs(gradMag(Y));
cnt = 64;
%GAVGx = filter2(window,GX,'valid')./cnt;
GAVGx_sqr = filter2(window,GX.*GX,'valid')./cnt - AVGx_sqr;
%GAVGy = filter2(window,GY,'valid')./cnt;
GAVGy_sqr = filter2(window,GY.*GY,'valid')./cnt - AVGy_sqr;
GAVGxy = filter2(window,GX.*GY,'valid')./cnt - AVGxy;

SSIM_MAP = ((2*AVGxy + c1).*(2*GAVGxy + c2))./((AVGx_sqr+AVGy_sqr+c1).*(GAVGx_sqr+GAVGy_sqr+c2));
MSSIM = mean(SSIM_MAP(:));

end

%Calculate mean using the integral image method
function M = integralMean(X)
I = integralImage(X);
winsz = 8;
w_end = size(I,2);
h_end = size(I,1);
wincnt = winsz*winsz;

%V1 = I(1:(h_end-winsz),1:(w_end-winsz));
%V2 = I(1:(h_end-winsz),(winsz+1):w_end);
%V3 = I((winsz+1):h_end,1:(w_end-winsz));
%V4 = I((winsz+1):h_end,(winsz+1):w_end);

M = I((winsz+1):h_end,(winsz+1):w_end)+I(1:(h_end-winsz),1:(w_end-winsz))-(I(1:(h_end-winsz),(winsz+1):w_end)+I((winsz+1):h_end,1:(w_end-winsz)));
M = M./wincnt;

end

%Generate a integral image
function I = integralImage(X)
 I = cumsum(cumsum(X,1),2);
end

% Calculate gradient
function G = gradMag(X)
    Rx = [1 0 ; 0 -1];
    Ry = [0 1 ; -1 0];
    
    Gx = abs(filter2(Rx,X,'valid'));
    Gy = abs(filter2(Ry,X,'valid'));
    
    G = max(Gy,Gx) + (1/4)*min(Gy,Gx);
end