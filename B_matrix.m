function b = B_matrix(projectors)
% Usage: b = B_matrix(projectors)
%
% Returns the b matrix for a given tomographically complete set of
% measurements.  This is defined by b_nu,mu =
% <psi_nu|gamma(mu)|psi_nu>

  dim_m = size(projectors, 2);
  dim_b = dim_m^2;
  tmp = zeros(dim_b);
  
  for i=1:dim_b
    for j=1:dim_b
      tmp(i,j) = projectors(i,:)*sigma_N(j, dim_m)*(projectors(i,:)');
    end
  end

  %why do the transpose???
  b = transpose(tmp);
