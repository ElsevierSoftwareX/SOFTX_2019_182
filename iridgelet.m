% X = IRIDGELET( RT, L )
% Computes inverse ridgelet transform as given in reference [1]. Suitable
% for denoising applications.
% 
% Input Parameters:
% RT = Ridgelet transform
% L = No. of scales in the ridgelet transform
% 
% Output Parameters:
% X = Reconstructed image
% 
% See also ridgelet
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

function x = iridgelet( rt, L, md )

if ~exist( 'md', 'var' )
  md = 0;
end

l = size( rt, 1 );
if md == 0
  yp = fft_iso_idwt( rt(1,:), L, md );
elseif md == 1
  yp = fft_iso_idwt( rt{1}, L, md );
end
for jj = 2:l
if md == 0
  yp(jj,:) = fft_iso_idwt( rt(jj,:), L, md );
elseif md == 1
  yp(jj,:) = fft_iso_idwt( rt{jj}, L, md );
end  
end
yr = rectopolar_2_cart( yp );
x = ifft2( yr );
[l1 l2] = size( x );
x = x(1:l1-1,1:l2-1);