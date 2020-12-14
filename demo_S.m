clc
clear
close all
addpath test_images solvers utilits

I = im2double(imread('barbara.png'));    %gray image
I = I(257:end,257:end);

% I = im2double(imread('kodim01.png'));   %rgb image
% I = I(201:456,129:128+256,:);

[n1,n2,n3] = size(I);
OPTK = 'S';

% Mask
K = floor(double(imread('mask3.bmp'))/255);
A = zeros(n1,n2,n3);
for i = 1:n3; A(:,:,i) = K; end
K = A;
clear A
x0   = K.*I;

opts.I  = I;
opts.verbose = 0;
opts.Tol   = 1e-3;
opts.MaxIt = 200;

opts.tau  = 0.001;  opts.mu = 0.015;
beta = 0.01;
opts.beta1 = beta;      opts.beta2 = beta;

% PPSM
alg_opts.gamma = 1.6;
alg_opts.s = 2.01;
alg_opts.r = 1;
[u,v,out] = CLRP_PPSM(x0,K,OPTK,opts,alg_opts);

figure;
subplot(2,2,1)
imshow(x0,[])
title('Observed')
subplot(2,2,2)
imshow(u,[])
title('Cartoon')
subplot(2,2,3)
imshow(v,[])
% imshow(rgb2gray(v),[])  % colorful texture
title('Texture')
subplot(2,2,4)
imshow(u+v,[])
title('Cartoon + Texture')
suptitle('PPSM')

% % EADM
% [u_EADM,v_EADM,out_EADM] = CLRP_EADM(x0,K,OPTK,opts);