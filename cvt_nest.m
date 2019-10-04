function [nstd nmean] = cvt_nest( nrow, ncol, J, L, md )
% [NSTD NMEAN] = CVT_NEST( NROW, NCOL, J, L, MD )
% Monte-Carlo estimation of the mean and standard deviation of the curvelet
% coefficients of white noise with size NROW X NCOL. The computed standard
% deviation may be used to determine threshold for image denoising, ref[1].
% 
% Input Parameters
% NROW, NCOL = Image size
% J = No. of levels in the wavelet pyramid in the start (no. of subbands in
%     the curvelet).
% L = A vector indicating the no. of scales in the local ridgelets; order =
%     fine to coarse. (e.g., when J = 4, L = [3 4 4 5]).
% 
% Output Parameters
% NSTD = Cell containing the standard deviation of the curvelet
% coefficients. To compute the standard deviation for each directional
% subband, use the function CELLMEAN.
% NMEAN = Cell containing the mean of the curvelet coefficients
% 
% References:
% [1] J.L. Starck, E.J. Candes, and D.L. Donoho, "The Curvelet Transform
% for Image Denoising," IEEE Trans. on Image Proc., Vol 11, No. 6, June
% 2002.
%
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: May 1, 2010
% Modified on: May 7, 2010

itrn = 30;

randn( 'state', rand(1)*1000 );
x = randn( nrow, ncol );
y = cvt( x, J, L, md );
nstd = cellsquare( y );
nmean = y;
disp(1);

for ii = 1:itrn-1
  randn( 'state', rand(1)*1000 );  
  x = randn( nrow, ncol );
  y = cvt( x, J, L, md );
  % sum of squares of the realizations
  nstd = celladd( nstd, cellsquare( y ) );
  % sum of the realizations
  nmean = celladd( nmean, y );
  disp(ii+1);
end

nstd = cellsubtract( nstd, celldivide( cellsquare( nmean ), itrn ) );
nstd = cellpower( celldivide( nstd, (itrn-1) ), 0.5 );
nmean = celldivide( nmean, itrn );