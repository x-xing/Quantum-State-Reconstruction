function out = reconstruct2q_li(counts, projectors, param)

if nargin == 2
    param(1) = 15;          %counting time
    param(2) = 190e-9;      %coincidence window
end

M = zeros(16,4,4);
rho = zeros(4,4);

B = B_matrix(projectors);
B_inv = inv(B);

n = counts(:,4) - counts(:,2).*counts(:,3)*param(2)/param(1)*4*1.04; 
for mu=1:1:16
    M(mu,:,:) = M_matrix(mu, projectors, B, B_inv);
end
    
for nu=1:1:16
    rho = rho + reshape(M(nu,:,:),4,4)*n(nu);
end
rho = rho/sum(n(1:4));
out = rho;