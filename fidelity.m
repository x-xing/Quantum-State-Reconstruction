function out = fidelity(rho1,rho2);
%function out = fidelity(rho1,rho2) calculates the fidelity between two
%matrix rho1 and rho2.

if size(rho1)== size(rho2)
    out=real(trace((rho1^0.5*rho2*rho1^0.5)^0.5));
else
    error('Dimensions of input matrices should be the same');
end