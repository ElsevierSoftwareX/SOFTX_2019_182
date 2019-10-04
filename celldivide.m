% Z = CELLDIVIDE(Z1,Z2)
% Subtracts the elements of the cells Z2 from Z1. If Z1 and Z2 contain subcells,
% this function calls itself recursively, untill it encounters a matrix or
% a number. The cells Z1 and Z2 should be number cells, and they should be
% of the same type and size; or Z2 can be a scalar.
% 
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: April 20, 2010
% Modified on: April 20, 2010

function z = celldivide(z1,z2)

[m n] = size(z1);

if ~all(size(z2) == [m n])
  if size(z2) ~= 1
    error('Cells should be of the same type and size or the second argument should be a scalar');
  else
    type = 1;
  end
else
  if all(size(z2) == 1)
    type = 1;
  else
  type = 0;
  end
end

z = cell(m,n);

if type == 0
  
  for jj = 1:m*n
    if iscell(z1{jj})
      z{jj} = celldivide(z1{jj},z2{jj});
    else
      z{jj} = (z1{jj}) ./ z2{jj};
    end
  end

else

  for jj = 1:m*n
    if iscell(z1{jj})
      z{jj} = celldivide(z1{jj},z2);
    else
      z{jj} = (z1{jj}) / z2;
    end
  end
  
end