function x = snr(sig, ref)
% Input:
%       sig         Modified image
%       ref         Reference image
%  
% Output:
%       x           SNR value
mse = mean((ref(:)-sig(:)).^2);
dv = var(ref(:),1);
x = 10*log10(dv/mse);
