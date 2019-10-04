% Y = CELLMEAN( X )
% Function to find the mean of the elements in a number cell. If a cell
% element is a cell itself, the function calls itself recursively, until a
% matrix or a number is met with.
% 
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: April 22, 2010
% Modified on: April 22, 2010

function y = cellmean( x )

[m n] = size( x );
y = cell( m, n );

for ii = 1:m*n
  if iscell( x{ii} )
    y{ii} = cellmean( x{ii} );
  else
    y{ii} = mean( x{ii}(:) );
  end
end