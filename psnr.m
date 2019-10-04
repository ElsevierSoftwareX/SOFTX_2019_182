% R = PSNR( X, XR, GRAYMAX )
% Returns the peak signal-to-noise ration (PSNR) between two images X and
% XR, where XR is the reconstruction of X using some method. GRAYMAX is the
% maximum gray value possible and is optional; the default value is 255.
% 
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: April 20, 2010
% Modified on: April 20, 2010.

function r = psnr( x, xr, graymax )

if ~exist( 'graymax', 'var' )
  graymax = 255;
end

r = 10*log10((graymax^2 * numel( x )) / norm( x-xr, 'fro' )^2);