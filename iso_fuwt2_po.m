function D = iso_fuwt2_po(x,J,hr)
% D = ISO_FUWT2_PO(X,J,Hr)
% 
% Performs 2D forward isometric undecimated discrete wavelet transform
% (alias stationary wavelet transform, maximum overlap DWT, A Trous
% algorithm etc.). In isometric UWT, the only one set of high
% frequency wavelet coeficients in a stage is computed by subtracting the
% current approximation coefficients from the previous approximation
% coefficients. Hence reconstruction is easy: just sum the wavelet
% coefficients from all levels and the approximation coefficients in the
% coarsest level.
% 
% Input Arguments:
% X  : input signal
% J  : No. of levels in the algorithm
% Hr : Scaling filter (lowpass in the synthesis side, 1D; 2D will be
%      computed by Kronecker tensor product)
% Output Arguments:
% D  : UWT coefficients; size(D) = [size(x), J+1]; 'J+1'st
% layer contains the Jth level approximation coefficients. In the third
% dimension, as the index increases, the coefficients goes from finer to
% coarse levels. See references.
% 
% See also iso_iuwt2_po
% 
% References: 
% [1] A Wavelet Tour of Signal Processing, Stephane Mallat, 2nd Ed.
% [2] J. Starck, J. Fadili and F. Murtagh, "The Undecimated Wavelet
%     Decomposition and its Reconstruction," IEEE Trans. on Imag. Proc., Vol
%     16, No.2, pp 297-309, Feb 2007.
% 
% Author: Sandeep Palakakl (sandeep.dion@mgial.com)
% Orgn: IIT Madras
% Created on: March 7, 2010

[ha] = on_wavelet_filters(hr);
ha = ha(:)';
ha = ha/sum(ha);
lf = length(hr);
delay = floor(lf/2);
% delay = lf-1;

[m n] = size(x);
K1 = ceil(log2(m));
K2 = ceil(log2(n));
if(m ~= 2^K1)
  disp('Warning: No. of rows of the signal not a power of 2; Padding zeros on top and down');
  ls = (2^K1-m);
  x = [zeros(floor(ls/2),n) x zeros(ceil(ls/2),n)];
end
if(n ~= 2^K2)
  disp('Warning: No. of columns of the signal not a power of 2; Padding zeros on front and back');
  ls = (2^K2-n);
  x = [zeros(m,floor(ls/2)) x zeros(m,ceil(ls/2))];
end

ha1 = ha'*ha;
hha1 = ha1;

A = x;
D = zeros(size(x));
for ii = 1:J
%   UWT: Lowpass
%   [mh nh] = size(hha1);
  A_new = fast_cconv2(cshift2(A,2^(ii-1)*(delay),2^(ii-1)*(delay),'u','l'),hha1);  
  D(:,:,ii) = A-A_new;
  A = A_new;
%   Preparation for the next iteration
  if(ii~=J)
    hha1 = upsample2(hha1,2);
  end
end
D(:,:,J+1) = A;