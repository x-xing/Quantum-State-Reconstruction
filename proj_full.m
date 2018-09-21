function projectors=proj_full;
% tt = 2.67;
% %tt = 2*pi/3;
% H=[cos(tt);-sin(tt)];
% V=[-sin(tt);-cos(tt)]; 
% 
H=[1;0];
V=[0;1];

R=(H-i*V)/sqrt(2);
L=(H+i*V)/sqrt(2);
D=(H+V)/sqrt(2);
A=(H-V)/sqrt(2);
projectors=zeros(36,4);

basis=[H,V,D,A,R,L];

for ind=1:36
    projectors(ind,:)=tensor_product(basis(:,fix((ind-1)/6)+1),basis(:,mod(ind-1,6)+1))';
end
