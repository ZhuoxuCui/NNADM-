    
function [sol,out,spsnr,relerr] = cadmimpulsivetv(A,y,opts,delta,x_true)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Min_{x,z,w} J^k_{lammada,rho1,rho2}(x,z,w)
% Input: 
%        y     the noise observation 
%    delta     noise level 
%    eigsA     eigen value of A
%¡¡¡¡   Aty¡¡   A'*y
%       CL     eigen value of AtA*rho1 + LtL*rho2
%  maxiter1     max number of iteration
%  maxiter2     K_max
%   x_true     true solution just need for psnr and snr
%   alpha      the parameter in f
%    beta      the parameter in g
% Output:
%     sol     solution of NNADM++
%    iter     number of iteration used in  NNADM++
%   spsnr     psnr histery
%  relerr     relative error histery
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ========================================================================
% Copyright (c), May, 2017
% Zhuo-Xu Cui  Dept. Math.,Whuhan Univiversity
% zhuoxucui@whu.edu.cn
% ========================================================================

% initialization for the main loop 
if nargin < 5; opts = []; end 
rho1 = opts.rho1; rho2 = opts.rho2; beta = opts.beta; alpha = opts.alpha;
eta = opts.eta; lammada=opts.lammada; gamma = opts.gamma;
maxiter2=opts.maxiter2;
maxiter1=opts.maxiter1;
tt = opts.tt;

sizeb = size(y);
eigsA = psf2otf(A,sizeb);
Aty = real(ifft2(conj(eigsA) .* fft2(y))); 
eigsDtD = abs(psf2otf([1,-1],sizeb)).^2 + abs(psf2otf([1;-1],sizeb)).^2;
eigsAtA = abs(eigsA).^2; 
x = real(ifft2(conj(eigsA) .* fft2(y)));

w = zeros(sizeb);
z1 = w;
z2 = z1;
D = @(U) Grad(U); % define Grad and Dive
Dt = @(X,Y) Dive(X,Y);
spsnr = [];       %  the vactor for store the psnr history 
relerr = [];
%% Main loop
iter = 1; % to count the seps used in the main loop 
iter1 = 1;
iter3=0;
ep = opts.ep;
for iter1=1:maxiter1
     CL = eta+rho1*eigsAtA + rho2*eigsDtD;
   for iter=1:maxiter2 
% ====================================================================================================      
    % ==================
    %     x-subprolem
    % ==================
    Atw = real(ifft2(conj(eigsA) .* fft2(w))); 
    gradphi2 = phiprim( real( ifft2(eigsA .* fft2(x)) )-y,beta,ep);
    [temp1,temp2] = D(x);
    normdx = sqrt(temp1.^2+temp2.^2);
    phitemp = phiprim(normdx,alpha,ep);
    right = rho1*(Aty+Atw)+rho2*Dt(z1,z2)+1/eta*x-real(ifft2(conj(eigsA) .* fft2(gradphi2)))-lammada*phitemp.*Dt(temp1,temp2);  
    x = fft2(right)./(CL);
    x = real(ifft2(x));
    psnrX = mpsnr(x,x_true);
    spsnr = [spsnr;psnrX];
    [D1x,D2x] = D(x);
    
    % ==================
    %    z-subprolem
    % ==================
    
    temp1 = D1x;
    temp2 = D2x;
    V = temp1.^2 + temp2.^2;
    V = sqrt(V);
    barV =  V;
    barV(V==0) = 1;
    V = max(V - lammada*alpha/rho2, 0)./barV;
    z1 = temp1.*V;
    z2 = temp2.*V;

    % ====================
    %     w-subprolem
    % ====================
    
    Axy = real( ifft2(eigsA .* fft2(x)) ) - y;
    V1 =  Axy.^2; 
    V1 = sqrt(V1);
    barV1 =  V1;
    barV1(V1==0) = 1;
    V1 = max(V1 - beta/rho1, 0)./barV1;
    w = Axy.*V1;
% ====================================================================================================    

    relerr1 = norm(x-x_true,'fro')/norm(x_true,'fro');
    relerr = [relerr;relerr1];

    % ==================================================
    %     modified Morozov's discrepancy principle
    % ==================================================
    
    residual = norm(Axy(:),1);
       if  residual <=  tt*delta;
          break;
       end
  end
     iter3 = iter3+iter;
     
    % =============================
    %     continuation technique
    % =============================
    
    rho1 = rho1/gamma;
    rho2 = rho2/gamma;
    lammada = lammada*gamma;
    beta = beta/gamma;
    alpha = alpha/gamma;
end
sol=x;
out = [];
out.iter3 = iter3;
out.iter1 = iter1;
out.rho1 = rho1;
out.rho2 = rho2;
out.lammada = lammada;
function y = phiprim(x, alpha_ep, ep)
         y = -ep*alpha_ep^2*(2+ep*alpha_ep*abs(x)) .*x ./ ((1+ep*alpha_ep*abs(x)).^2);
end

end