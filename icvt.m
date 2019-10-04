% X = ICVT( CV, J, L )
% Computes the inverse curvelet (version 1) of an image, as given in [1].
% 
% Input Parameters:
% CV = Curvelet coefficients.
% J = No. of levels in the wavelet pyramid in the start (no. of subbands in
%     the curvelet).
% L = A vector indicating the no. of scales in the local ridgelets; order =
%     fine to coarse.
% Output Parameters:
% X = Input image.
%
% See also cvt local_ridgelet local_iridgelet
%
% References:
% [1] J.L. Starck, E.J. Candes, and D.L. Donoho, "The Curvelet Transform
% for Image Denoising," IEEE Trans. on Image Proc., Vol 11, No. 6, June
% 2002.
%
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: May 1, 2010
% Modified on: May 7, 2010

function x = icvt( cv, J, L, md )

if ~exist( 'md', 'var' )
  md = 0;
end

% Min local window size of the local ridgelet tfm
Bmin = 16;
% Inverse local ridgelet
B = Bmin;
D(:,:,J+1) = cv{1};
for jj = 2:J+1
  if numel( L ) > 1
    lev = L(jj-1);
  else lev = L;
  end
  D(:,:,jj-1) = local_iridgelet( cv{J+3-jj}, B, lev, md );
  if rem(jj-1,2) == 1
    B = 2*B;
  end
end
% Inverse IUWT
x = iso_iuwt2_po(D);