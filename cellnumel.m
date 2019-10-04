% Z = CELLNUMEL(Z)
% Computes the total number of primary elements in a cell.
% 
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: Sep 21, 2010
% Modified on: Sep 21, 2010

function N = cellnumel(z)

  N = 0;
if iscell(z)
  [m n] = size(z);
  for jj = 1:m*n
    if iscell(z{jj})
      N = cellnumel(z{jj});
    else
      N = N + numel(z{jj});
    end
  end
else
  N = N + numel(z);
end