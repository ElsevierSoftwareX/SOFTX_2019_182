function A = fdct_wrapping_dispcoef_expand(u,v,B)
  A = zeros(u,v);
  [p,q] = size(B);
  A(1:p,1:q) = B;