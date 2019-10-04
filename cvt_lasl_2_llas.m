% Y = CVT_LASL_2_LLAS( X )
% Curvelet coefficients are saved in the order:
% Level (IUWT), Location (Block), Angle, Scale (Ridgelet).
% This is referred to as llas.
% Another arrangement for denoising algorithms is:
% Level (IUWT), Angle, Scale (Ridgelet), Location (Block).
% This is referred to as lasl
% This function changes the arrangment from lasl to llas.
% 
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: June 2010

function y = cvt_lasl_2_llas( x )

nlev = length( x );
y{1} = x{1};
for ii = 2:nlev
  nangle = length( x{ii} );
  nscale = length( x{ii}{1} );
  nloc = size( x{ii}{1}{1}, 1 );
  y{ii} = cell( 1, nloc );
  for jj = 1:nloc
    y{ii}{jj} = cell( 1, nangle );
    for kk = 1:nangle
      for mm = 1:nscale
        y{ii}{jj}{kk}{mm} = x{ii}{kk}{mm}(jj,:);
      end
    end
    y{ii}{jj} = y{ii}{jj}';
  end
end