MLHFlag=0; %Maximum likelihood not implemented yet.

proj4q=zeros(256,16);
M = zeros(256,16,16);
global n;
n = zeros(256,1);
rho = zeros(16,16);
counting_timef = 10;
coin_win = 160e-9;
errload = 0;

%%%% Construct the 4-qubit projectors from 2-qubit projectors. This is a
%%%% bit time consuming, but it only need to be done once.
proj2q=proj_james;
proj2q_p = proj_path('ul');
for ind1=1:16
    for ind2=1:16
            proj4q(ind1*16-16+ind2,:)=tensor_product(proj2q(ind1,:),proj2q_p(ind2,:));
        end
end
save('proj4q_pre_cal.mat','proj4q');

%%%% Calculate the B Matrix. Only need to be done once.
B = B_matrix(proj4q);
save('B_pre_cal.mat','B');
B_inv = inv(B);
save('B_inv_pre_cal.mat','B_inv');

%%%% Calculate the M Matrix. Only need to be done once.
load('B_pre_cal.mat');
load('B_inv_pre_cal.mat');
for mu=1:1:256
    M(mu,:,:) = M_matrix(mu, proj4q, B, B_inv);
end
save('M_pre_cal.mat','M');

load('proj4q_pre_cal.mat');
load('M_pre_cal.mat');

%%%% Construct the theoretical counts, n
%%1. Theoretical rho for 1/sqrt(2)*(|HHHH>+ i |VVVV>)
% rho_t([1 16],[1 16])=[0.5 0.5i;-0.5i 0.5]; 
% %2. Theoretical rho for loop cluster
qinit = [1 0]';
Hadamad =   [1 1; 1 -1]/sqrt(2);
q_plus = qo(Hadamad*qinit);
graph_ini = tensor(tensor(q_plus,q_plus),tensor(q_plus,q_plus));
Q1 = qgate(4,cphase,[1,2]);
Q2 = qgate(4,cphase,[2,4]);
Q3 = qgate(4,cphase,[1,4]);
Q4 = qgate(4,cphase,[1,3]);
lpGraph = Q4*Q3*Q2*Q1*graph_ini;
rho_t = double(lpGraph)*double(lpGraph)';
for ind1=1:256
        n_t(ind1)=100*proj4q(ind1,:)*rho_t*proj4q(ind1,:)';  %assuming there are 100 4-fold coincidences.
end
n=n_t;

p1 = pwd;
folder = [p1(1:strfind(p1, 'My Dropbox')-1) 'My Dropbox\Data\Tomo\2011-06-19\16Basis\'];
dirListing = dir([folder '2011-06-19_afterInterfero*']);
for d = 1:length(dirListing)
    fileName = fullfile(folder,dirListing(d).name); 
    counts_1file = dlmread(fileName,'\t');
    s1f=counts_1file(:,6);
    s2f=counts_1file(:,7);
    s3f=counts_1file(:,9);
    s4f=counts_1file(:,13);
    c13f=counts_1file(:,10)-s1f.*s3f.*coin_win.*4./counting_timef*0.92; 
    c14f=counts_1file(:,14)-s1f.*s4f.*coin_win.*4./counting_timef*0.92; 
    c23f=counts_1file(:,11)-s2f.*s3f.*coin_win.*4./counting_timef*0.92; 
    c24f=counts_1file(:,15)-s2f.*s4f.*coin_win.*4./counting_timef*0.92; 
    c13f(find(c13f<0))=0; c14f(find(c14f<0))=0; c23f(find(c23f<0))=0; c24f(find(c24f<0))=0;
    if size(findstr(dirListing(d).name,'HH'),2) == 1
        %pol_proj(1,2,3,4) = ch23, ch24, ch14, ch13
        %n(pol_proj*16-16+(1:16))
        n(1*16-16+(1:16)') = c23f;
        n(2*16-16+(1:16)') = c24f;
        n(3*16-16+(1:16)') = c14f;
        n(4*16-16+(1:16)') = c13f;     
        errload = errload + 1;
    end        
    if size(findstr(dirListing(d).name,'HD'),2) == 1
        %pol_proj(12,13) = ch23, ch13
        n(12*16-16+(1:16)') = c23f;
        n(13*16-16+(1:16)') = c13f;
        errload = errload + 1;
    end  
    if size(findstr(dirListing(d).name,'HR'),2) == 1
        %pol_proj(14,15) = ch14, ch24
        n(14*16-16+(1:16)') = c14f;
        n(15*16-16+(1:16)') = c24f;
        errload = errload + 1;
    end  
    if size(findstr(dirListing(d).name,'DH'),2) == 1
        %pol_proj(7,8) = ch24, ch23
        n(7*16-16+(1:16)') = c24f;
        n(8*16-16+(1:16)') = c23f;
        errload = errload + 1;
    end  
    if size(findstr(dirListing(d).name,'DD'),2) == 1
        %pol_proj(10) = ch23
        n(10*16-16+(1:16)') = c23f;
        errload = errload + 1;
    end  
    if size(findstr(dirListing(d).name,'DR'),2) == 1
        %pol_proj(9) = ch23
        n(9*16-16+(1:16)') = c23f;
        errload = errload + 1;
    end  
    if size(findstr(dirListing(d).name,'RH'),2) == 1
        %pol_proj(5,6) = ch23, ch24
        n(5*16-16+(1:16)') = c23f;
        n(6*16-16+(1:16)') = c24f;
        errload = errload + 1;
    end  
    if size(findstr(dirListing(d).name,'RD'),2) == 1
        %pol_proj(11) = ch23
        n(11*16-16+(1:16)') = c23f;
        errload = errload + 1;
    end  
    if size(findstr(dirListing(d).name,'RR'),2) == 1
        %pol_proj(16) = ch24
        n(16*16-16+(1:16)') = c24f;
        errload = errload + 1;
    end      
end

if errload < 9
    error('Data loading not successful. Some basis missing.');
else if errload >9
        error('Data loading not successful. Multiple files found in the same basis.');
     end        
end
ndisp = reshape(n,16,16);


for nu=1:1:256
    rho = rho + reshape(M(nu,:,:),16,16)*n(nu);
end

rho = rho/sum(n([1:4,17:20,33:36,49:52])); %Renomalization, in H, V basis. HHHH, HHHV, etc.

f1=figure;
subplot(1,2,1);
bar3(real(rho));
title('real part');

subplot(1,2,2);
bar3(imag(rho));
title('imaginary part');

d1=eig(rho);
purity1 = sum(diag(rho*rho));

f2=figure; bar3(abs(rho))
title('Amplitude (Linear Inversion)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Maximum likelyhood estimation
if MLHFlag
    t=FindInitialT(rho);
    fhandle=@fun_MLH;
    
    [t,fval]=fminsearch(fhandle,t,optimset('MaxIter',600*length(t),'MaxFunEvals',600*length(t)));
    rho_mlh=fun_rho(t);
    
    d2=eig(rho_mlh);
    purity2 = sum(diag(rho_mlh*rho_mlh));
    
    figure; bar3(abs(rho_mlh))
   % set(gca,'XTickLabel',{'HH','HV','VH','VV'},'YTickLabel',{'HH','HV','VH','VV'});
    title('Amplitude (MLH Estimation)');
    disp([,'Eig_Linear','   Eig_MLH']);
    disp(real([d1,d2]));
    disp(['Purity_Linear','  Purity_MLH']);
    disp(real([purity1,purity2]));
else
    disp('Eig_Linear');
    disp(real(d1));
    disp('Purity_Linear');
    disp(real(purity1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




