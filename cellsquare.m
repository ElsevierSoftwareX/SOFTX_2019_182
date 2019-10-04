% Z = CELLSQUARE(Z)
% Squares the elements of the cell Z. If Z contain subcells, this function
% calls itself recursively, untill it encounters a matrix or a number. The
% cell should contain only numbers.
% 
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: April 20, 2010
% Modified on: April 20, 2010

function z = cellsquare(z)

[m n] = size(z);
for jj = 1:m*n
  if iscell(z{jj})
    z{jj} = cellsquare(z{jj});
  else
    z{jj} = (z{jj}).^2;
  end
end