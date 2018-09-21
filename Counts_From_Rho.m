clear;
phi=pi;
rho_des=1/2*[0 0 0 0;0 1 exp(i*phi) 0; 0 exp(-i*phi) 1 0; 0 0 0 0];

n=zeros(1,16);
msg='';
proj1=proj_james;
for ind=1:16
    n(ind)=sum(diag(rho_des*proj1(ind,:)'*proj1(ind,:)))*1000;
    msg=sprintf('%s;%3.2f',msg,n(ind));
end
disp(msg);




    