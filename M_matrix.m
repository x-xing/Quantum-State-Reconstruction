function m = M_matrix(mu, projectors, B, B_inv)
% Usage: m = M_matrix(mu, projectors, B, B_inv)
%
% Returns the mu'th M matrix used for quantum state tomography

  dim_m = size(projectors, 2);
  dim_b = dim_m^2;

  if(nargin < 3)
    B=chop(B_matrix(projectors));
  end
  if(nargin < 4)
    B_inv = inv(b);
  end
  
  tmp = zeros(size(projectors, 2));
  for j=1:dim_b
    tmp = tmp+B_inv(mu,j) * sigma_N(j, dim_m);
  end
  m = tmp;
