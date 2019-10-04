% XP = CART_2_RECTOPOLAR( XR )
% Returns the center slices of signal by assuming that the input is a
% 2D-DFT, and slicing is done by mapping the input to 2D-Fourier transform
% plane, using fftshift for 2D. The output is again mapped to 1D-DFT axis,
% using ifftshift in 1D. Slicing is done as given in [1]. This is used for
% computing radon transform, ridgelet transform and curvelet transform.
% 
% Input Parameters:
% XR = Data in rectangular coordinate; square size NxN; even or odd.
%
% Output Parameters:
% XP = Data in polar coordinate, size: 4*floor(N/2) x N
%
% See also rectopolar_2_cart ridgelet iridgelet
%
% References:
% [1] J.L. Starck, E.J. Candes, and D.L. Donoho, "The Curvelet Transform
% for Image Denoising," IEEE Trans. on Image Proc., Vol 11, No. 6, June
% 2002.
%
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: April 30, 2010
% Modified on: May 7, 2010

function xp = cart_2_rectopolar(xr)

[m n] = size(xr);
if m ~= n
  error( 'Input array to rectopolar should be square' )
end

xr = fftshift(xr);

oe = 0;
if rem( n, 2) == 0
  xr(n+1,:) = 0;
  xr(:,n+1) = 0;
  oe = 1;
  n = n+1;
end

flh = floor(n/2);
xp = zeros(4*flh,n);
for l = 1:2*flh
    ml = (flh-l)/flh;
    x = -flh:flh;
    y = round(ml*x)+flh+1;
    for k = 1:n
      xp(l,k) = xr(k,y(k));
      xp(l+n-1,k) = xr(n+1-y(k),k);      
    end
end

xp = cshift1(xp,2*flh ,'d');
% xp = cshift1(xp,flh+2 ,'d');

if oe == 1
  xp = xp(:,1:n-1);
end

xp = ifftshift(xp,2);