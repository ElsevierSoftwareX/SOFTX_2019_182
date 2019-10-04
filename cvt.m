% CV = CVT( X, J, L )
% Computes the curvelet (version 1) of an image, as given in [1].
% 
% Input Parameters:
% X = Input image.
% J = No. of levels in the wavelet pyramid in the start (no. of subbands in
%     the curvelet).
% L = A vector indicating the no. of scales in the local ridgelets; order =
%     fine to coarse. (e.g., when J = 4, L = [3 4 4 5]).
% 
% Output Parameters:
% CV = Curvelet coefficients.
%
% See also icvt local_ridgelet local_iridgelet
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

function cv = cvt( x, J, L, md )

if ~exist( 'md', 'var' )
  md = 0;
end

% Min local window size of the local ridgelet tfm.
Bmin = 16;
% IUWT of the image
hr = (1/16)*[1 4 6 4 1];
D = iso_fuwt2_po( x, J, hr ); % fine to coarse
% Local ridgelet
cv = cell(1,J+1);
B = Bmin;
lev = L;
for jj = 1:J
  if numel( L ) > 1
    lev = L(jj);
  end
  cv{J+2-jj} = local_ridgelet( D(:,:,jj), B, lev, md );
  if rem(jj,2) == 1
    B = 2*B;
  end
end
cv{1} = D(:,:,J+1); % cv coarse to fine