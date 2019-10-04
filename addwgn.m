function [y wn] = addwgn(x,snr)
% function [Y WN] = ADDWGN(X,SNR)
% X: 1D input  signal.
% SNR: SNR required in dB.
% Y: the output signal.
% WN: Added noise (optional argument).
%
% Adds white Gaussian noise (AWGN) of zero mean to 1D signal x so that the
% signal-to-noise ration is snr in dB. The zero-mean and unit-variance WGN
% is computed by randn, and then it is scaled to maintain the required SNR
% given by 'snr' in the function call. 
% 
% Author: Sandeep P, IITM
% Created on: Oct 12, 2009
% Modified on: Oct 26, 2009

xe = sum(x(:).^2); % Energy of the signal
% Noise
w = randn(size(x));
w = w-mean(w(:));
we = sum(w(:).^2); % The energy of the noise
sig = sqrt(xe./(10^(snr/10)*we)); % scaling coefficient to get the desired snr.
w = sig*w;
y = x + w;
if(nargout == 2)
    wn = w;
end

% EOF