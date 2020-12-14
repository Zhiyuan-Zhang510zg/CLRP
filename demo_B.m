clc
clear
close all
addpath test_images solvers utilits

I = im2double(imread('barbara.png'));   %gray image
I = I(1:256,257:end);

% I = im2double(imread('kodim03.png'));   %rgb image
% I = I(1:64*8,129:64*10,:);

[n1,n2,n3] = size(I);
opts.MaxIt = 200;
OPTK = 'B';
opts.I = I;
Blur_Radius = 5;
K = fspecial('disk',Blur_Radius);
x0 = imfilter(I,K,'circular');

opts.Tol = 1e-2;
opts.verbose = 0;

opts.tau = 1e-5;  opts.mu = 1e-4;
beta = 1e-3;
opts.beta1 = beta;      opts.beta2= beta;

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
title('Texture')
subplot(2,2,4)
imshow(u+v,[])
title('Cartoon + Texture')
suptitle('PPSM')

% % EADM
% [u_EADM,v_EADM,out_EADM] = CLRP_EADM(x0,K,OPTK,opts);


