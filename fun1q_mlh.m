function outL = fun1q_mlh(t,n_prob,proj)
Rho = t*t';

dim_p = size(proj, 1);
tmp=zeros(1,dim_p);

for ind=1:dim_p
    tmp(ind) = proj(ind,:)*Rho*(proj(ind,:)');
end

outL=sum((tmp-n_prob).^2/2./double(tmp));