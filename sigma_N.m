function matrix=sigma_N(j, N)
% Usage: matrix=sigma_N(j, N)
%
% Returns the j'th generalized generator for SU(N).  sigma_N(N^2) is
% always eye(N), the rest of the generators are normalized so that
% trace(sigma(i) * sigma(j)) = N*delta(i,j)

if(j<1 || j>N^2)
  error('sigma_n: j out of range for SU(N)');
end

m = fix((j-1)/N)+1;
n = mod((j-1), N)+1;
tmp1 = zeros(N,1);
tmp2 = tmp1;
tmp1(m) = 1;
tmp2(n) = 1;

if(m<n)
  matrix = (tmp1*tmp2' + tmp2*tmp1')*sqrt(N/2);
elseif(m>n)
  matrix = i*(tmp1*tmp2' - tmp2*tmp1')*sqrt(N/2);
elseif(m < N)
  matrix = -sqrt(N/(m^2+m)) * diag(1:N <= m);
  matrix(m+1, m+1) = m*sqrt(N/(m^2+m));
else % m==n==N
  matrix = eye(N);
end
