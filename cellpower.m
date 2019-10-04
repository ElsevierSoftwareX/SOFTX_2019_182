% Z = CELLPOWER(Z,Q)
% Finds the power of the elements of the cell Z to a real number Q. If Z
% contain subcells, this function calls itself recursively, untill it
% encounters a matrix or a number. The cell should contain only numbers.
% 
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: April 20, 2010
% Modified on: April 20, 2010

function z = cellpower(z,q)

[m n] = size(z);
for jj = 1:m*n
  if iscell(z{jj})
    z{jj} = cellpower(z{jj},q);
  else
    z{jj} = (z{jj}).^q;
  end
end