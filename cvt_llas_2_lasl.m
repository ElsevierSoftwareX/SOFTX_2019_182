% Y = CVT_LLAS_2_LASL( X )
% Curvelet coefficients are saved in the order:
% Level (IUWT), Location (Block), Angle, Scale (Ridgelet).
% This is referred to as llas.
% Another arrangement for denoising algorithms is:
% Level (IUWT), Angle, Scale (Ridgelet), Location (Block).
% This is referred to as lasl
% This function changes the arrangment from llas to lasl.
% 
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: June 2010

function y = cvt_llas_2_lasl( x )

nlev = length( x );
y = cell( 1, nlev );
y{1} = x{1};
for ii = 2:nlev
  nloc = length( x{ii} );
  nangle = length( x{ii}{1} );
  nscale = length( x{ii}{1}{1} );
  y{ii} = cell( 1, nangle );
  for jj = 1:nangle
    y{ii}{jj} = cell( 1, nscale );
    for kk = 1:nscale
      y{ii}{jj}{kk}{1} = x{ii}{1}{jj}{kk};
      for mm = 2:nloc
        y{ii}{jj}{kk}{mm} = x{ii}{mm}{jj}{kk};
      end
      y{ii}{jj}{kk} = cell2mat(y{ii}{jj}{kk}(:));
    end
  end
end