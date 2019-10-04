% Y = CSHIFT2(X,M,N,UD,RL)
% Performs a 2D circular shift on a matrix x on the
% columns and/or rows.
% UD = 'u' : up-shift of columns;
%      'd' : down-shift of columns;
% RL = 'l' : left-shift of rows;
%      'r' : right-shift of rows;
% 
% Y = CSHIFT2(X,M,RLUD)
% Performs a 2D circular shift on a matrix x on the
% columns or rows, according to RLUD.
% RLUD = 'r' : right-shift of rows;
%        'l' : left-shift of rows;
%        'u' : up-shift of columns;
%        'd' : down-shift of columns;
% 
% See also cshift1
% 
% Author: Sandeep P (email: sandeep.dion@gmail.com)
% Orgn: IIT Madras
% Created on: Feb 23, 2010
% Modified on: March 2, 2010

function y = cshift2(x,m,n,ud,rl)

if(nargin == 3)
  if(n == 'l' | n == 'r')
    rl = n;
    n = m;
    m = 0;
    ud = 'u';
  elseif(n == 'u' | n == 'd')
    ud = n;
    n = 0;
    rl = 'l';
  end
end

rlud1 = ud;
rlud2 = rl;
[p q] = size(x);
if(rlud1 == 'd')
  m = -mod(m,p);
elseif(rlud1 == 'u') m = mod(m,p);
  else disp('Error in rlud1')
    return
end
if(rlud2 == 'r')
  n = -mod(n,q);
elseif(rlud2 == 'l') n = mod(n,q);
  else disp('Error in rlud2')
    return
end
y = zeros(p,q);
% left-up shift
y = x(mod(m:p+m-1,p)+1,mod(n:q+n-1,q)+1);