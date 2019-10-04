% Y = UPSAMPLE2(X,M,N)
% Performs upsampling of 2D functions. X should be a matrix (not a vector).
% M and N are the upsampling factors along colmns (downward direction) and
% rows (horizontal direction), respectively.
% 
% Y = UPSAMPLE2(X,M)
% Same as above; N = M;
% 
% Author: Sandeep Palakkal (email: sandeep.dion@gmail.com)
% Created on: Feb 23, 2010
% Modified on: Feb 23, 2010

function y = upsample2(x,m,n)

if(nargin == 2)
  n = m;
end
y = upsample(x',n);
y = upsample(y',m);