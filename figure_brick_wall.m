function figure_brick_wall
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
downsampleFactor=4;
I = imread('brick_wall.tiff');
   if downsampleFactor>1
        I = imresize(I,1/downsampleFactor);
   end
I = double(I);
x_true = I/max(max(I));
framekd = 1;
Nlev = 1;
opts = [];
opts.ep = 1;
opts.rho1 = 1000;
opts.rho2 = 10;
opts.lammada = 0.02;
opts.alpha = 10;
opts.beta = 1;
opts.eta  = 1;
opts.tt = 1.001;
opts.gamma = 0.8;
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
opts.maxiter2 = 70;
thtype = 0;
%%  precompute the constant in the main loop
disp('--------------NNADM++ is running------------')
tic,
[sol,out,spsnr,relerr] = cadmimpulsiveframe(A,y,opts,delta,framekd,Nlev,x_true,thtype);
t1 = toc
%%
relerror=norm(sol-x_true,'fro')/norm(x_true,'fro'); 
iter = length(spsnr)
%% Plot result
figure(3); imshow(sol);
title(sprintf('NNADM++, PSNR %4.2fdB, CPU %4.2fs',mpsnr(sol,x_true),t1),'fontsize',40);
fprintf('PSNR(y) %4.2fdB, PSNR(Recovered) %4.2fdB,',mpsnr(y,x_true),mpsnr(sol,x_true))
fprintf(' Iteration %d\n\n',iter)
figure(4)
semilogy(1:iter, spsnr,'k-','LineWidth', 2);
title('PSNR')
xlabel('iteration'); ylabel('PSNR');