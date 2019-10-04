% LRT = LOCAL_RIDGELET( X, B, L, md )
% Computes local ridgelet as given in [1]. The ridgelet function is also
% implemented as given in [1].
% 
% Input Parameters:
% X = Input image
% B = Block size
% L = No. of scales in the ridgelet
% md = if 1, o/p is cell; if 0, o/p is a vector.
%
% Output Parameters:
% LRT = Local ridgelet transform coefficients.
%
% See also local_iridgelet ridgelet iridgelet
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

function lrt = local_ridgelet( x, B, L, md )

if ~exist( 'md', 'var' )
  md = 0;
end

[m n] = size( x );
if m~=n
  error( 'Input must be a square array' );
end

bl = 2*n/B;
ll = bl*bl;
lrt = cell(1,ll);
for ii = 1:bl-1
  ri = (ii-1)*(B/2)+1:(ii+1)*(B/2);
  for jj = 1:bl-1
    ci = (jj-1)*(B/2)+1:(jj+1)*(B/2);
    y = x(ri,ci);
    poi = ll+1-(ii-1)*bl-jj;
    lrt{1,poi} = ridgelet( y, L, md );
  end
  ci = [n-(B/2)+1:n 1:B/2];
  y = x(ri,ci);
  poi = ll+1-ii*bl;
  lrt{1,poi} = ridgelet( y, L, md );
end

ri = [n-(B/2)+1:n 1:B/2];
for jj = 1:bl-1
  ci = (jj-1)*(B/2)+1:(jj+1)*(B/2);
  y = x(ri,ci);
  poi = ll+1-(bl-1)*bl-jj;
  lrt{1,poi} = ridgelet( y, L, md );
end
ci = [n-(B/2)+1:n 1:B/2];
y = x(ri,ci);
poi = ll+1-bl*bl;
lrt{1,poi} = ridgelet( y, L, md );