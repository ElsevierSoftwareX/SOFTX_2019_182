% RT = RIDGELET( X, L, md )
% Computes ridgelet transform as given in reference [1]. Suitable for
% denoising applications.
% 
% Input Parameters:
% X = Input image of square size
% L = No. of scales in the ridgelet transform
% md = if 1, o/p is cell; if 0, o/p is a vector.
% 
% Output Parameters:
% RT = Ridgelet transform of X
% 
% See also iridgelet
%          fft_iso_dwt         fft_iso_idwt
%          cart_2_rectopolar   rectopolar_2_cart
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

function rt = ridgelet( x, L, md )

if ~exist( 'md', 'var' )
  md = 0;
end

[l1 l2] = size( x );
if rem( l1,2 ) == 0  && rem( l2,2 ) == 0
  pt1 = l1 + 1; pt2 = l2 + 1;
end

yr = fft2( x, pt1, pt2 );
yp = cart_2_rectopolar( yr );
l = size( yp, 1 );
if md == 0
  rt = fft_iso_dwt( yp(1,:), L, md );
elseif md == 1
  rt = cell(l,1);
  rt{1} = fft_iso_dwt( yp(1,:), L, md );
end
for jj = 2:l
  if md == 0
    rt(jj,:) = fft_iso_dwt( yp(jj,:), L, md );
  elseif md == 1
    rt{jj} = fft_iso_dwt( yp(jj,:), L, md );
  end
end