% Y = CELLSIGN( X )
% Applies the signum function on a cell X. If X contains subcells, then the
% function works recursively. The output is a cell of the same structure as
% X, and conatins the sign of the elements of X.
% 
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: April 20, 2010
% Modified on: May 20, 2010

function y = cellsign( x )

[m n] = size( x );

y = cell( m, n );

for jj = 1:m*n
  if iscell( x{jj} )
    y{jj} = cellsign( x{jj} );
  else
      y{jj} = sign( x{jj} );
  end
end