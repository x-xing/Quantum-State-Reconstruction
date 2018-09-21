function projectors=proj_2q_pol
% tt = 2.67;
% %tt = 2*pi/3;
% H=[cos(tt);-sin(tt)];
% V=[-sin(tt);-cos(tt)]; 

H=[1;0];
V=[0;1];
%%%%%%%%%%%%%%%%%%%%%%%%%
%due to negative hwp & qwp angle, R<->L, D<->A
 L=(H-i*V)/sqrt(2);     
R=(H+i*V)/sqrt(2);
A=(H+V)/sqrt(2);
D=(H-V)/sqrt(2);
%%%%%%%%%%%%%%%%%%%%%%%%%
projectors=zeros(16,4);

%operator '  : the complex conjugated, transpose array
%operator .' : the non-complex conjugated, transpose array
%We should use the complexed transpose here???

projectors(1,:) =tensor_product(H,H)';
projectors(2,:) =tensor_product(H,V)';
projectors(3,:) =tensor_product(V,V)';  
projectors(4,:) =tensor_product(V,H)';
projectors(5,:) =tensor_product(R,H)';
projectors(6,:) =tensor_product(R,V)';
projectors(7,:) =tensor_product(D,V)';
projectors(8,:) =tensor_product(D,H)';
projectors(9,:) =tensor_product(D,R)';
projectors(10,:)=tensor_product(D,D)';
projectors(11,:)=tensor_product(R,D)';
projectors(12,:)=tensor_product(H,D)';
projectors(13,:)=tensor_product(V,D)';
projectors(14,:)=tensor_product(V,L)';
projectors(15,:)=tensor_product(H,L)';
projectors(16,:)=tensor_product(R,L)';