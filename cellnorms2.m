% Y = CELLNORMS2( X )
% Computes the l2-norm of the individual number subcells in a number cell
% X, by considering each subcell as a vector. The cell X should be a number
% cell.
%
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: April 20, 2010
% Modified on: Feb 15, 2011

function y = cellnorms2( x )

if iscell( x )
  [m n] = size( x );
  y = cell( m, n );
  for jj = 1:n
    for kk = 1:m
      y{(kk-1)*m+jj} = cellnorms2( x{(kk-1)*m+jj} );
    end
  end

else
  y = sqrt(sum(abs(x(:)).^2));
end