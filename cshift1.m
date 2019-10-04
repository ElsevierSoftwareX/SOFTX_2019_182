% Y = CSHIFT1(X,M,RLUD)
% Performs a 1D circular shift of length m on a vector of matrix x on the
% columns or rows, depending up on rlud.
% rlud = 'l' : left-shift of rows;
%        'r' : right-shift of rows;
%        'u' : up-shift of columns;
%        'd' : down-shift of columns;
% 
% See also cshift2
% 
% Author: Sandeep P (email: sandeep.dion@gmail.com)
% Created on: Feb 23, 2010
% Modified on: Feb 23, 2010

function y = cshift1(x,m,rlud);

[p q] = size(x);
y = zeros(p,q);
if(rlud == 'l')
  len = q;
  m = mod(m,len);
  y(:,1:len-m) = x(:,m+1:len);
  y(:,len-m+1:len) = x(:,1:m);
elseif(rlud == 'r')
  len = q;
  m = mod(m,len);
  y(:,1:m) = x(:,len-m+1:len);  
  y(:,m+1:len) = x(:,1:len-m);
elseif(rlud == 'u')
  len = p;
  m = mod(m,len);  
  y(1:len-m,:) = x(m+1:len,:);
  y(len-m+1:len,:) = x(1:m,:);
elseif(rlud == 'd')
  len = p;
  m = mod(m,len);
  y(m+1:len,:) = x(1:len-m,:);
  y(1:m,:) = x(len-m+1:len,:);
end