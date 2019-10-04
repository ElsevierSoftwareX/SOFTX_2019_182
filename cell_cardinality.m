% Y = CELL_CARDINALITY( X )
% Computes the total number of elements in a cell X. If X contains
% subcells, this function calls itself recursively, untill it encounters a matrix or
% a number. The cell X should be a number cell.
% 
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: April 20, 2010
% Modified on: Dec 21, 2010

function y = cell_cardinality( x )

[m n] = size( x );

y = 0;

for jj = 1:m*n
  if iscell( x{jj} )
    y = cell_cardinality( x{jj} ) + y;
  else
    y = numel( x{jj} ) + y;
  end
end