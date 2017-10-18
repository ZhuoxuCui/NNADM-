function figure_cameraman
% ========================================================================
% Copyright (c), May, 2017
% Zhuo-Xu Cui
% zhuoxucui@whu.edu.cn 
% ========================================================================
clc
clear all;
close all;
addpath('solvers/');
addpath('solvers/coresolvers/');
addpath('solvers/utilities/');
path(path,genpath(pwd));
addpath(genpath(fileparts(mfilename('fullpath'))));
warning off 
%%
downsampleFactor=1;
I = imread('camera256.bmp');
   if downsampleFactor>1
        I = imresize(I,1/downsampleFactor);
   end
I = double(I);
x_true = I/max(max(I));
opts = [];
opts.ep = 1;
opts.rho1 = 800;
opts.rho2 = 5;
opts.lammada = 0.02;
opts.alpha = 10;
opts.beta = 1;
opts.eta  = 1;
opts.tt = 1.001;
opts.gamma = 0.85;
A = fspecial('gaussian',[20 20], 30); 
%A = fspecial('motion',50,90);
%% add blur and noise 
noisetype = 'saltpepper'; % 'saltpepper' or 'gaussian'
d_per     = 0.2;          % percentage of noise
randn('state',0)
y1 = imfilter(x_true,A,'circular','conv');  
if strcmp(noisetype,'saltpepper')
   y = imnoise(y1, 'salt & pepper',d_per);
end
noise = y1-y;
delta = norm(noise(:),1) % the real noise level 
rate = delta/norm(y);
figure(1);
imshow(x_true);
title('Original','fontsize',40);
figure(2);
imshow(y);
title(sprintf('Observed, PSNR %4.2fdB',mpsnr(y,x_true)),'fontsize',40);
%% set parameter 
opts.maxiter1 = 30;
opts.maxiter2 = 40;
%%  precompute the constant in the main loop
disp('--------------NNADM++ is running------------')
tic,
[sol,out,spsnr,relerr] = cadmimpulsivetv(A,y,opts,delta,x_true);
t1 = toc

%%
relerror=norm(sol-x_true,'fro')/norm(x_true,'fro'); 
iter = length(spsnr)
%% Plot result
figure(3); imshow(sol);
title(sprintf('NNADMM++, PSNR %4.2fdB, CPU %4.2fs',mpsnr(sol,x_true),t1),'fontsize',40);
fprintf('PSNR(y) %4.2fdB, PSNR(Recovered) %4.2fdB,',mpsnr(y,x_true),mpsnr(sol,x_true))
fprintf(' Iteration %d\n\n',iter)