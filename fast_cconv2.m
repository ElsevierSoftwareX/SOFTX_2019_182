% Y = FAST_CCONV2(X,H,OPT)
% Performs 2D circular convolution. This implementation is fast since it
% uses FFT2 and IFFT2.
% 
% Input Parameters:
% X  = Input Image.
% H  = 2D filter.
% OPT = 'null': No effect, (default option)
%       'nc_shift': treats H as a noncausal filter with origin at
%                   floor(size(H)/2); the output is shifted by the above
%                   amount in up-and-left direction to correct the origin.
% 
% Output Parameters:
% Y  = Circular convolution of X and H. No. of rows of Y is the maximum of
%      the number of rows of X and that of H. Simlilarly the no. of columns
%      of Y.
% 
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Date: April 13, 2010
% Modified on: April 13, 2010

function y = fast_cconv2(x,h,opt)
  
if(~exist('opt','var')) 
    opt = 'null';
end

[mx nx] = size(x);
[mh nh] = size(h);
m = max(mx,mh);
n = max(nx,nh);
hf = fft2(h,m,n);
xf = fft2(x,m,n);
y = real(ifft2(xf.*hf));
switch opt
  case 'null'
%     Do nothing
  case 'nc_shift'
  msh = floor(mh/2);
  nsh = floor(nh/2);
  y = cshift2(y, msh, nsh, 'u', 'l');
  otherwise
    disp('Unknown option for "opt" in fast_cconv2');
end