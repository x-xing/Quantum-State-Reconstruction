function out = reconstruct1q_li(n_prob, projectors)

M = zeros(4,2,2);
rho = zeros(2,2);

B = B_matrix(projectors);
B_inv = inv(B);

for mu=1:1:4
    M(mu,:,:) = M_matrix(mu, projectors, B, B_inv);
end
    
for nu=1:1:4
    rho = rho + reshape(M(nu,:,:),2,2)*n_prob(nu);
end
rho = rho/sum(n_prob(1:2));
out = rho;