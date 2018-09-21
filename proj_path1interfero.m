function projectors=proj_path1Interfero
% Projectors for the path qubit. Syntax: proj_path(a)
% Use a string input 'a' to specify the interferometer used. 
% a = 'uu', 'ul', 'lu', or 'll'.
% Note for the path qubit, the transmission of the tomography PBS is logical 1.

% tt = 2.67;
% %tt = 2*pi/3;
% H=[cos(tt);-sin(tt)];
% V=[-sin(tt);-cos(tt)]; 
% 
% if a == 'uu'
%     phi1=0; phi2=0;
% else if a =='ul'
%         phi1=0; phi2=-0.9342;
%     else if a == 'lu'
%             phi1=0.0986; phi2=0;
%         else if a=='ll'
%                 phi1 =  0.0986;phi2 = -0.9342;
%             else 
%                 error('Wrong input paramater for proj_path(a). a= ''uu'', ''ul'', ''lu'', or ''ll''.');
%             end
%         end
%     end
% end

phi1 = 0;
phi2 = 0;

V1=[1;0];
H1=[0;exp(i*phi1)];
R1=(H1-i*V1)/sqrt(2);
L1=(H1+i*V1)/sqrt(2);
D1=(H1+V1)/sqrt(2);
A1=(H1-V1)/sqrt(2);

V2=[1;0];
H2=[0;exp(i*phi2)];
R2=(H2-i*V2)/sqrt(2);
L2=(H2+i*V2)/sqrt(2);
D2=(H2+V2)/sqrt(2);
A2=(H2-V2)/sqrt(2);

projectors=zeros(16,4);

%operator '  : the complex conjugated, transpose array
%operator .' : the non-complex conjugated, transpose array
%We should use the complexed transpose here???


projectors(1,:) =tensor_product(H1,H2)';
projectors(2,:) =tensor_product(H1,V2)';
projectors(3,:) =tensor_product(V1,V2)';  
projectors(4,:) =tensor_product(V1,H2)';
projectors(5,:) =tensor_product(R1,H2)';
projectors(6,:) =tensor_product(R1,V2)';
projectors(7,:) =tensor_product(D1,V2)';
projectors(8,:) =tensor_product(D1,H2)';
projectors(9,:) =tensor_product(D1,R2)';
projectors(10,:)=tensor_product(D1,D2)';
projectors(11,:)=tensor_product(R1,D2)';
projectors(12,:)=tensor_product(H1,D2)';
projectors(13,:)=tensor_product(V1,D2)';
projectors(14,:)=tensor_product(V1,L2)';
projectors(15,:)=tensor_product(H1,L2)';
projectors(16,:)=tensor_product(R1,L2)';