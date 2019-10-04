% Xr = ISO_IUWT2_PO(D)
% 
% Performs 2D inverse isometric undecimated discrete wavelet transform
% (alias stationary wavelet transform, maximum overlap DWT, A Trous
% algorithm etc.). In isometric UWT, the only one set of high
% frequency wavelet coeficients in a stage is computed by subtracting the
% current approximation coefficients from the previous approximation
% coefficients. Hence reconstruction is easy: just sum the wavelet
% coefficients from all levels and the approximation coefficients in the
% coarsest level.
% 
% Input Arguments:
% D  : UWT coefficients; size(D) = [size(x), J+1]; 'J+1'st
% layer contains the Jth level approximation coefficients. In the third
% dimension, as the index increases, the coefficients goes from finer to
% coarse levels. See references.
% Output Arguments:
% X  : Reconstructed input signal
% 
% See also iso_fuwt2_po
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

function  xr = iso_iuwt2_po(D)

% Reconstruction
xr = sum(D,3);