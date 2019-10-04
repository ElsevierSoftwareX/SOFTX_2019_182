% MSS = CELLCOMPARE( Y, CTH )
% Compares the elements of the cells Y and CTH. CTH can be a scalar. If Y
% contains subcells, then comparison is performed recursively. If both are
% cells, then they should have the same structure. The output MSS is a cell
% with the same structure as Y, and the elements are 1, if Y is greater
% than CTH at that position, and 0, otherwise.
%
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: April 22, 2010
% Modified on: April 22, 2010

function MSS = cellcompare( y, cth )

if iscell(y)
  [m n] = size( y );
  if iscell(cth)
    if ~all( size( cth ) == [m n] )
      error('Cells should be of the same type and size or the second argument should be a scalar');
    else
      type = 0;
    end
  else
    if numel(cth) == 1
      type = 1;
    else
      error('The second argument should be either a cell of the same size as the first one or a scalar');
    end
  end
end

MSS = cell( m, n );

if type == 0

  for jj = 1:m*n
    if iscell( y{jj} )
      MSS{jj} = cellcompare( y{jj}, cth{jj} );
    else
      MSS{jj} = (y{jj}) > cth{jj};
    end
  end

else

  for jj = 1:m*n
    if iscell( y{jj} )
      MSS{jj} = cellcompare( y{jj},cth );
    else
      MSS{jj} = (y{jj}) > cth;
    end
  end

end