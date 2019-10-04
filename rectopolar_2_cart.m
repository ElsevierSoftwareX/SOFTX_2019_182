% XR = RECTOPOLAR_2_CART( XP )
% Rectopolar to Cartesian coordinate conversion. Inverts the effect of
% the function cart_2_rectopolar
% Input Parameters:
% XP = Data in polar coordinate.
%
% Output Parameters:
% XR = Data in rectangular coordinate; square size NxN; even or odd.
%
% See also cart_2_rectopolar ridgelet iridgelet
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

function xr = rectopolar_2_cart( xp )

[m n] = size(xp);
if m ~= 4*floor(n/2)
  error( 'Size mismatch on input' )
end

xp = fftshift(xp,2);

oe = 0;
if rem( n, 2) == 0
  xp(:,n+1) = 0;
  oe = 1;
  n = n+1;
end

flh = floor(n/2);
xp = cshift1(xp,2*flh ,'u'); 
% xp = cshift1(xp,flh+2 ,'u');

xr = zeros( n, n );
sidx = zeros( n, n );
for l = 1:2*flh
  ml = (flh-l)/flh;
  x = -flh:flh;    
  y = round(ml*x)+flh+1;
  for k = 1:n
    xr(k,y(k)) = xr(k,y(k)) + xp(l,k);
    sidx(k,y(k)) = sidx(k,y(k)) + 1;    
    xr(n+1-y(k),k) = xr(n+1-y(k),k) + xp(l+n-1,k);
    sidx(n+1-y(k),k) = sidx(n+1-y(k),k) + 1;
  end
end

if oe == 1
  xr = xr(1:n-1,1:n-1);
  sidx = sidx(1:n-1,1:n-1);
end

xr = xr ./ sidx;
xr = ifftshift(xr);