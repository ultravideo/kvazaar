% Test different step-wise downscaling schemas

mex -g downScaler.c

%% Load image and define other param

img = imread('saturn.png');%imread('peppers.png');
img = rgb2ycbcr(img);
s =(size(img(:,:,1)));

trgts = (0.9:-0.1:0.1);
trgts = uint32(repmat(trgts,2,1).*repmat(s',1,size(trgts,2)));

%% Baseline (orig -> target_)
% Loop over targets
loop = 1;
br = cell(1,size(trgts,2));
bg = br;
bb = br;
bt = 0; %total time
avgt = 20; %how many times the scaling is done per trgt to get a better time
for trgt = trgts
    t = 0;
    for i = 1:avgt
        tic;
        [br{loop},bg{loop},bb{loop}] = downScaler( img(:,:,1), trgt', img(:,:,2), trgt', img(:,:,3), trgt');
        t = toc;
    end
    bt = bt + t/avgt;
    loop = loop+1;
end
disp(['Time taken for in total for baseline: ' num2str(bt)]);

%% step-wise (orig->target1->target2->...)
% Loop over targets
loop = 2;
sr = {img(:,:,1)};
sg = {img(:,:,2)};
sb = {img(:,:,3)};
st = 0; %total time
avgt = 20; %how many times the scaling is done per trgt to get a better time
for trgt = trgts
    t = 0;
    for i = 1:avgt
        tic;
        [sr{loop},sg{loop},sb{loop}] = downScaler( sr{loop-1}, trgt', sg{loop-1}, trgt', sb{loop-1}, trgt');
        t = toc;
    end
    st = st + t/avgt;
    loop = loop+1;
end
disp(['Time taken for in total for step-wise: ' num2str(st)]);

%% Check if there is error between the methods
disp(['Time diff ratio: ' num2str(st/bt*100) '%']);

num_pics = size(trgts,2);
error = zeros(3,num_pics);
er_map = cell(3,num_pics);
for l = 1:size(trgts,2)
    [error(1,l),er_map{1,l}] = ssim(sr{l+1},br{l});
    [error(2,l),er_map{2,l}] = ssim(sg{l+1},bg{l});
    [error(3,l),er_map{3,l}] = ssim(sb{l+1},bb{l});
end

disp('Errors: ');
disp(error);

% Display error maps
for l = 1:num_pics
    subplot(3,num_pics,l);
    imshow(er_map{1,l});
    title(['R trgt' num2str(l)]);
    subplot(3,num_pics,l+num_pics);
    imshow(er_map{2,l});
    title(['G trgt' num2str(l)]);
    subplot(3,num_pics,l+2*num_pics);
    imshow(er_map{3,l});
    title(['B trgt' num2str(l)]);
end