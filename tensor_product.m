function rv=tensor_product(varargin)
% Usage: AB = tensor_product(A, B, ...)
%
% Calculates the tensor product of two matricies/vectors A and B.
% This is a matrix whose elements are all possible products of an
% element from A and an element from B in the form of 
%
% [a_11*B a_12*B a_13*B ...
%  a_21*B a_22*B a_23*B ...
%  ... ];
%
% This is useful on state vector, density matricies, and to tile a
% matrix in a repetative form.
%
% More than two input arguments can be provided, in which case they
% are evaluated as:
%
% AxBxC = (AxB)xC

rv = 1;
for j = 1:nargin
  next = varargin{j};
  [n11, n12] = size(rv);
  [n21, n22] = size(next);
  
  ind11=(1:n11)';
  ind12=(1:n12)';
  ind21=(1:n21)';
  ind22=(1:n22)';

  ind11=ind11(:,ones(1,n21));
  ind12=ind12(:,ones(1,n22));
  ind21=ind21(:,ones(1,n11));
  ind22=ind22(:,ones(1,n12));

  rv = rv(ind11',ind12') .* next(ind21,ind22);
end
