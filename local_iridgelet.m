% X = LOCAL_IRIDGELET( LRT, B, L )
% Computes inverse local ridgelet as given in [1]. The ridgelet function is
% also implemented as given in [1].
% 
% Input Parameters:
% LRT = Local ridgelet transform coefficients.
% B = Block size
% L = No. of scales in the ridgelet
%
% Output Parameters:
% X = Reconstructed image
%
% See also local_ridgelet ridgelet iridgelet
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

function x = local_iridgelet( lrt, B, L, md )

if ~exist( 'md', 'var' )
  md = 0;
end

ll = size( lrt, 2 );
bl = sqrt( ll );
n = bl*B/2;
win = (cos(pi*(-B/2:B/2-1)/(B))).^2;
win = win'*win;
x = zeros( n, n );
for ii = 1:bl-1
  ri = (ii-1)*(B/2)+1:(ii+1)*(B/2);
  for jj = 1:bl-1
    ci = (jj-1)*(B/2)+1:(jj+1)*(B/2);    
    poi = ll+1-(ii-1)*bl-jj;    
    y = iridgelet( lrt{1,poi}, L, md );
    x(ri,ci) = x(ri,ci) + y.*win;
  end
  ci = [n-(B/2)+1:n 1:B/2];
  poi = ll+1-ii*bl;
  y = iridgelet( lrt{1,poi}, L, md );
  x(ri,ci) = x(ri,ci) + y.*win;
end

ri = [n-(B/2)+1:n 1:B/2];
for jj = 1:bl-1
  ci = (jj-1)*(B/2)+1:(jj+1)*(B/2);
  poi = ll+1-(bl-1)*bl-jj;
  y = iridgelet( lrt{1,poi}, L, md );
  x(ri,ci) = x(ri,ci) + y.*win;
end
ci = [n-(B/2)+1:n 1:B/2];
poi = ll+1-bl*bl;
y = iridgelet( lrt{1,poi}, L, md );
x(ri,ci) = x(ri,ci) + y.*win;