% X = FFT_ISO_IDWT( WC, L )
% Inverse isometric wavelet transform using B3 spline function as scaling filter
% frequcny response. The scaling filter is thus compact in the frequency
% domain. The transformation is done in the frequency domain using FFT. The
% function returns a wavelet pyramid that has twice as many samples as in
% the input sequency.
% 
% Input Parameter:
% WC = Wavelet coefficients of length twice that of the signal, arranged
%      from coarse to fine. Finest level is the last N samples, if N is the
%      signal length. The next coarser level is the final N/2 samples in
%      the remaining part, and so on. If L is full (L = log2(N)), then the
%      coarsest component is the first sample and the second sample is
%      zero.
% L = No. of wavelet levels.
% 
% Output Parameter:
% X = DFT of the reconstructed signal.
% 
% See also fft_iso_dwt.m
% 
% References:
% [1] J.L. Starck et al., Image Processing and Data Analysis: The
%     Multiscale Approach, Cambridge, U.K., 1998.
% [2] J.L. Starck et al., "Image reconstruction by the wavelet transform
%     applied to aperture synthesis," Astron. Astrophys., vol. 283, pp.
%     349-360, 1994.
%
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: April 29, 2010.
% Modified on: May 7, 2010.

function xf = fft_iso_idwt( wc, L, md )

if exist( 'md', 'var' )
  if md == 1
    wc = cell2mat( wc );
  end
end


wcl = length( wc );
if 2^(floor(log2( wcl ))) == wcl
  oe = 0;
  xln = wcl/2;
else
  oe = 1;
  xln = floor((wcl-L-1)/2)+1;
end
  
xlev = floor(log2( xln ));
if 2^xlev+oe ~= xln
  error( 'Unsuitable input length' )
end
x = wc(1:2^(xlev-L+1)+oe);
wc = wc(2^(xlev-L+1)+oe+1:end);
xf = fft(x);
for jj = 1:L
  wf = fft( wc(1:2^(xlev-L+jj)+oe) );
  wc = wc(2^(xlev-L+jj)+oe+1:end);
  xf = xf + wf;
  if jj ~= L
    xf = xf(1:floor(length(xf)/2)+1);
    xf = [xf zeros(1,length(xf)-1)]; %#ok<AGROW>
    xf = [xf conj(fliplr(xf(2:end-1+oe)))];  
  end
end