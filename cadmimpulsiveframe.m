    
function [sol,out,spsnr,relerr] = cadmimpulsiveframe(A,y,opts,delta,framekd,Nlev,x_true,thtype)
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
%   beta       the parameter in g
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
[n1,n2] = size(y);
if framekd == 0
    COFF=[0 1;1 1];
elseif framekd == 1
    COFF = [0 1 1;1 1 1;1 1 1];
else
    COFF = [1 1 1 1 1;1 1 1 1 1;1 1 1 1 1;1 1 1 1 1;1 1 1 1 1];
end
FrameD = @(U) FraDecML(U,framekd,Nlev); % define multilevel decompositon operator 
FrameS = @(U) FraRecML(U,framekd,Nlev); % define multilevel sythesis  operator

if nargin < 5; opts = []; end 
rho1 = opts.rho1; rho2 = opts.rho2; beta = opts.beta; alpha = opts.alpha;
eta = opts.eta; lammada=opts.lammada; gamma = opts.gamma;
maxiter2=opts.maxiter2;
maxiter1=opts.maxiter1;
tt = opts.tt;

sizeb = size(y);
eigsA = psf2otf(A,sizeb);
Aty = real(ifft2(conj(eigsA) .* fft2(y))); 
eigsLtL = 1;
eigsAtA = abs(eigsA).^2; 
x = real(ifft2(conj(eigsA) .* fft2(y)));
w = zeros(sizeb);
z = FrameD(x);
spsnr = [];       %  the vactor for store the psnr history 
relerr = [];
%% Main loop
iter = 1; % to count the seps used in the main loop 
iter1 = 1;
iter3=0;
ep = opts.ep;
for iter1=1:maxiter1
     CL = eta+rho1*eigsAtA + rho2*eigsLtL;
   for iter=1:maxiter2 
% ====================================================================================================      
    % ==================
    %     x-subprolem
    % ==================
    Atw = real(ifft2(conj(eigsA) .* fft2(w))); 
    gradphi2 = phiprim( real( ifft2(eigsA .* fft2(x)) )-y,beta,ep);%grad g_1
    temp = FrameD(x);
    normlx = coefnorm( temp,Nlev,COFF,n1,n2);% ||Lx||_2
    phitemp = phiprim(normlx,alpha,ep);
    right = rho1*(Aty+Atw)+rho2*FrameS(z)+1/eta*x-real(ifft2(conj(eigsA) .* fft2(gradphi2)))-lammada*phitemp.*FrameS(temp);  
    x = fft2(right)./(CL);
    x = real(ifft2(x));
    psnrX = mpsnr(x,x_true);
    spsnr = [spsnr;psnrX];   
    lx = FrameD(x);
    
    % ==================
    %    z-subprolem
    % ==================
    
    if thtype  == 1
      z = istropshrink(lx,COFF,Nlev,lammada*alpha/rho2,rho2);
    else
      z = anistropicshrink(lx,lammada*alpha/rho2,rho2); 
    end

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
       if  residual <=  tt*delta
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
out.beta = beta;
function y = phiprim(x, alpha_ep, ep)
         y = -ep*alpha_ep^2*(2+ep*alpha_ep*abs(x)) .*x ./ ((1+ep*alpha_ep*abs(x)).^2);
end

function  val  = coefnorm(M,Nlev,COFF,n1,n2)
        [mC,nC] = size(COFF);
        S = zeros(n1,n2);
        for ki = 1:Nlev
            for ji = 1:mC
              for jj = 1:nC
                 S = S + (M{ki}{ji,jj}).^2;
              end
            end
        end 
        val =  sqrt(S);
end

function R = anistropicshrink(A,th,lp)
%=========================================================
%  Execute the adaptive shrinkage of Framelet coefficients
%  stored in array A with threshhold tau. The low
%  "frequencies" coefficients are keep without any
%  changes.
%========================================================
L = length(A);
R = A;
for l = 1:L
   [NH,NW] = size(A{l});
   for nh = 1:NH
      for nw = 1:NW
          %temp = A{l}{nh,nw}*(lp*th);% since (lp*th) is near 1 ;
          temp = A{l}{nh,nw};  
          if ((nw == 1) && (nh == 1))
             R{l}{nh,nw} = temp;
          else
             R{l}{nh,nw} = sign(temp).*max(0,abs(temp)-th/(2^(l-1)));
          end   
      end
   end
end
end

function W = istropshrink(tempW,COFF,Nlev,th,lp)
%===================================================
%  Execute the adaptive group shrinkage of Framelet 
%  coefficients  stored in array tempW with threshhold
%  tau. The low  "frequencies" coefficients are keep 
%   without any  changes.
%====================================================
       W = tempW;
       [mC,nC] = size(COFF);
        for ki = 1:Nlev
            normtempW = 0;
            for ji = 1:mC
              for jj = 1:nC
                 normtempW = (COFF(ji,jj)*tempW{ki}{ji,jj}).^2 + normtempW;
              end
            end
            normtempW = normtempW.^0.5;
            normtempW(normtempW == 0) = 1;
            normtempW = max(normtempW - th/(2^(ki-1)), 0)./normtempW;
            for ji=1:mC
              for jj=1:nC
                  temp = tempW{ki}{ji,jj};
                  if COFF(ji,jj)==0
%                  W{ki}{ji,jj} = temp*(th*lp); % since (lp*th) is near 1 
                   W{ki}{ji,jj} = temp;
                  else
                   W{ki}{ji,jj} =  temp.*normtempW;
                  end
              end
            end
        end
end
end