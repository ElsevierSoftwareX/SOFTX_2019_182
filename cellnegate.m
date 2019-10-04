% Y = CELLABS( X )
% Negates the elements of the cell X. If X contains subcells, this function
% calls itself recursively, untill it encounters a matrix or a number. The
% cell X should be a number cell.
% 
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: April 20, 2010
% Modified on: May 21, 2010

function y = cellnegate( x )

[m n] = size( x );

y = cell( m, n );

for jj = 1:m*n
  if iscell( x{jj} )
    y{jj} = cellnegate( x{jj} );
  else
      y{jj} = ~( x{jj} );
  end
end