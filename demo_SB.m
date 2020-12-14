clc
clear
close all
addpath test_images solvers utilits

% Barbara
Current_Image_Name = 'barbara';
I = im2double(imread('barbara.png'));
I = I(1:256,257:end);

% % Barbarargb_RGB
% Current_Image_Name = 'barbarargb';
% I = im2double(imread('barbarargb.png'));
% I = I(11:10+512,221:220+512,:);

% % Kodim_wall
% Current_Image_Name = 'kodim01';
% I = im2double(imread('kodim01.png'));
% I = I(201:456,351:606,:);

% % Notty
% Current_Image_Name = 'notty';
% I = double(imread('notty.jpg'))/255;
% I = I(1:256,1:256,:);


[n1,n2,n3] = size(I);
opts.MaxIt = 200;
opts.I    = I;
opts.Tol = 1e-2;
opts.verbose = 0;

% Blurring kernel
if n1==512
    Blur_Radius = 3;
elseif n1==256
    Blur_Radius = 5;
end
K    = fspecial('disk',Blur_Radius);

% Mask 1
if n1==512 && n3==3
    Missing_rate = 0.3;
elseif n3 == 3
    Missing_rate = 0.4;
elseif n3 == 1
    Missing_rate = 0.3;
end
dr = randperm(n1*n2);
mesnum =round(n1*n2*(1-Missing_rate));
OM = dr(1:mesnum);
S = zeros(n1,n2);
S(OM) = 1;

% Mask 2
if strcmp(Current_Image_Name,'kodim01')
    rng(2020)
    dtx= 128; dty = 128; S1 = rand(dtx,dty)>0.2; S2 = ones(n1/dtx,n2/dty);
    S  = kron(S1,S2);
    Missing_rate = 1-sum(S(:))/numel(S);
end

% Mask 3
if strcmp(Current_Image_Name,'notty')
    S = floor(double(rgb2gray(imread('block_text.bmp')))/255);
    S = S(25:32*8+24,9:32*8+8);
    Missing_rate = 1-sum(S(:))/numel(S);
end

A =zeros(size(I));
for i = 1:n3
    A(:,:,i) = S;
end
S = A;
clear A
x0   = imfilter(I,K,'circular');
x0 = x0.*S;

opts.tau  = 1e-4;  opts.mu = 2e-4;
beta = 5e-2;
opts.beta1 = beta;      opts.beta2= beta;

% PPSM
alg_opts.gamma = 1.6;
alg_opts.s = 2.01;
alg_opts.r = 1;
[u,v,out] = CLRP_PPSM_SB(x0,K,S,opts,alg_opts);

figure;
subplot(2,2,1)
imshow(x0,[])
title('Observed')
subplot(2,2,2)
imshow(u,[])
title('Cartoon')
subplot(2,2,3)
if n3 == 3
    imshow(rgb2gray(v),[])
else
    imshow(v,[])
end
title('Texture')
subplot(2,2,4)
imshow(u+v,[])
title('Cartoon + Texture')
suptitle('PPSM')

% % EADM
% [u_EADM,v_EADM,out_EADM] = CLRP_EADM_SB(x0,K,S,opts);
