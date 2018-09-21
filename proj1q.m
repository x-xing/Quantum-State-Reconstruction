function out = proj1q(a)
H=[1;0];
V=[0;1];
b=sqrt(1-a^2);

D=b*H+a*V;
A=a*H-b*V;
L=b*H+i*a*V;
R=a*H-i*b*V;
proj=zeros(4,2);

proj(1,:) = H';
proj(2,:) = V';  
proj(3,:) = D';
proj(4,:) = A';
proj(5,:) = L';
proj(6,:) = R';
out = proj;