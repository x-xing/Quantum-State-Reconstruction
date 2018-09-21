function val = gamma(mu, nqubits)
% Usage: val = gamma(mu, nqubits)
%
% Returns the mu'th gamma matrix for quantum state tomography.
% These are the appropriately normalized generators of SU(2)^N
  
%  ind = multiloop_index(mu, 4*ones(1,nqubits))-1;
%  val = sigma(ind)/2;
  val = sigma_N(mu, 2^nqubits)/sqrt(2);
