% WC = FFT_ISO_DWT( X, L )
% Forward isometric wavelet transform using B3 spline function as scaling
% filter frequcny response. The scaling filter is thus compact in the
% frequency domain. The transformation is done in the frequency domain
% using FFT. The function returns a wavelet pyramid that has twice as many
% samples as in the input sequency.
% 
% Input Parameter:
% X = DFT of the input signal.
% L = No. of wavelet levels.
% 
% Output Parameter:
% WC = Wavelet coefficients of length twice that of the signal, arranged
%      from coarse to fine. Finest level is the last N samples, if N is the
%      signal length. The next coarser level is the final N/2 samples in
%      the remaining part, and so on. If L is full (L = log2(N)), then the
%      coarsest component is the first sample and the second sample is
%      zero.
% 
% See also fft_iso_idwt.m
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

function wc = fft_iso_dwt( xf, L, md )

if ~exist( 'md', 'var' )
  md = 0;
end

xln = length( xf );
xlev = floor(log2( xln ));
oe = 0;
if mod(xln,2^xlev) == 1
  oe = 1;
elseif mod(xln,2^xlev) ~= 0
  error( 'Unsuitable signal length in the input' )
end

if L > xlev
  error( 'No. of levels must in the wavelet transform must be smaller than the dyadic level of the signal')
end

% hf = abs( wave_spline_fil( xln ) );
hf = wave_b3spline_fil( xln );
hfj = hf;
cf = xf .* hfj;
wf = xf - cf;
if md == 0
  wc = ifft( wf );
elseif md == 1
  wc = cell( 1, L+1 );
  wc{L+1} = ifft( wf );
end

if(L ~= 1)
  xf = cf(1:floor(xln/4)+1);
  xf = [xf conj(fliplr(xf(2:end-1+oe)))];
  hfj = ifftshift(downsample(fftshift(hfj),2));
else
  xf = cf;
end

for jj = 2:L
  cf = xf .* hfj;
  wf = xf - cf;
  w = ifft( wf );
  if md == 0
    wc = [w wc];
  elseif md == 1
    wc{L+2-jj} = w;
  end
  
  if jj ~= L
    xf = cf(1:floor(xln/2^(jj+1))+1);
    xf = [xf conj(fliplr(xf(2:end-1+oe)))];
    hfj = ifftshift(downsample(fftshift(hfj),2));
  else
    xf = cf;
  end
end
if md == 0
  wc = [ifft(xf) wc];
elseif md == 1
  wc{1} = ifft( xf );
end