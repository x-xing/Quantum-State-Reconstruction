function outL=fun_MLH_4q(t)
global n;
Rho = fun_rho(t);

N=sum(n(1:4));
if size(n,2)==16
    projectors=proj_james;
else if size(n,2)==36
        projectors=proj_full;
    else
        msgbox('Error in input count data. Neither 16 or 36!');
    end
end

dim_p = size(projectors, 1);
tmp=zeros(size(n,1),size(n,2));

for i=1:dim_p
    tmp(i) = projectors(i,:)*Rho*(projectors(i,:)');
end

outL=sum((N*tmp-n).^2/2./(N*tmp));
