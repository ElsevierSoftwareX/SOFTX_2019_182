% HF = WAVE_B3SPLINE_FIL(PT)
% Returns the scaling filter frequency response for wavelet derived from
% the B3-spline as scaling filter frequency response. This scaling function
% is compact in the frequency domain.
% 
% Input Parameters:
% PT  = No. of fft points in the output HF; should be an integer above or
% equal to 16, preferably a power of 2 (E.g., 16, 32, 64, 128 ... ).
% 
% Output Parameters:
% H = Frequency response of the scaling filter. It is real, and hence
%     zero-phase.
% 
% References:
% [1] J.L. Starck et al., Image Processing and Data Analysis: The
%     Multiscale Approach, Cambridge, U.K., 1998.
% [2] J.L. Starck et al., "Image reconstruction by the wavelet transform
%     applied to aperture synthesis," Astron. Astrophys., vol. 283, pp.
%     349-360, 1994.
%
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: April 30, 2010.
% Modified on: April 30, 2010.

function z = wave_b3spline_fil(pt)

if pt < 16
  error( 'wave_b3spline_fil works only for input arguments above or equal to 16' );
end

ln = floor(pt / 8);
x1 = ones( 1, 2*ln );
y1 = conv( x1,x1 );
z1 = conv( y1,y1 );
x2 = ones( 1,ln );
y2 = conv( x2,x2 );
z2 = conv( y2,y2 );
M = ceil( length( z2 )/2 );
N = ceil( length( z1 )/2 );
za = z1 / max( z1 );
zb = z2 / max( z2 );
z = zb ./ za(N-M+1:N+M-1);
L = ceil( length( z )/2 );
z = [zeros( 1, ceil(pt/2)-L+1 ) z zeros( 1, floor(pt/2)-L )];
z = fftshift( z );