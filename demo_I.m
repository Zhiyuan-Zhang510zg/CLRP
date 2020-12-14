clc
clear
close all
addpath test_images solvers utilits

%% Mixed Image Test
opts.tau =0.05; opts.mu = 1; beta = 1;      % Mixed Image

% % Boy + texture4
cart = im2double(imread('boy.tif'));
text = im2double(imread('texture4.tif'));
I = 0.7*cart+0.3*text;

% % TomAndJerry + wool
% cart = im2double(imread('TomAndJerry.png'));
% text = im2double(imread('wool.png'));
% I = 0.7*cart+0.3*text;

% % Synthetic image
% I  =Synthesis_circles_fun(64,3,0.00);

% %% Natural Image test
% opts.tau  = 0.016; opts.mu = 0.2; beta = 2;% Natural Image

% % Mandril
% I = im2double(imread('mandril_gray.tif'));
% 
% 
% % Barbara
% 
% I = im2double(imread('barbara.png'));
% I = I(1:256,257:end);


[n1,n2,n3] = size(I);
opts.MaxIt = 100;
OPTK = 'I';
K    = 1;
opts.I     = I;
x0   = I;
S = 1;

opts.beta1 = beta; opts.beta2 = beta;
opts.verbose = 0;
opts.Tol  = 1e-2;

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