function out=makephysical(rho)
%function out=makephysical(rho)make sure rho is positive definite. 
% sloppy way of make positive by deleting the negative eigenvalues.
        [v d]=eig(rho);
        d(find(d<0))=0;
        d=d/trace(d);
        rho = v*d/v;
        out=tril(rho,-1)+ tril(rho,-1)'...
            +real(diag(diag(rho)))+eye(size(rho))*1e-6;
        
