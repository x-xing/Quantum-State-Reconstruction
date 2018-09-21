clear all;
MLHFlag=1; %Maximum likelihood not implemented yet.

proj4q=zeros(256,16);
M = zeros(256,16,16);
n = zeros(256,1);
rho = zeros(16,16);
rho_graph = zeros(16,16);
rho_tr = zeros(16,16);

counting_time = 15;
coin_win = 190e-9;
errload = 0;

%%%% Construct the 4-qubit projectors from 2-qubit projectors. This is a
%%%% bit time consuming, but it only need to be done once.
% proj1=proj_2q_pol;
% proj2 = proj_2q_pol;    %necessary to apply a sigmax operation later.
% for ind1=1:16
%     for ind2=1:16
%             proj4q(ind1*16-16+ind2,:)=tensor_product(proj1(ind1,:),proj2(ind2,:));
%         end
% end
% save('proj4q_pre_cal.mat','proj4q');
% 
% %%%% Calculate the B Matrix. Only need to be done once.
% B = B_matrix(proj4q);
% save('B_pre_cal.mat','B');
% B_inv = inv(B);
% save('B_inv_pre_cal.mat','B_inv');
% 
% %%%% Calculate the M Matrix. Only need to be done once.
% load('B_pre_cal.mat');
% load('B_inv_pre_cal.mat');
% for mu=1:1:256
%     M(mu,:,:) = M_matrix(mu, proj4q, B, B_inv);
% end
% save('M_pre_cal.mat','M');
disp('Loading pre-calculated projectors...');
load('proj4q_pre_cal.mat');
load('M_pre_cal.mat');

p1=pwd;
cd([p1(1:strfind(p1, 'Dropbox')-1) 'Dropbox\MatlabCode\utility']);
qopath;
cd(p1);
%%%% Construct the theoretical counts, n
%%1. Theoretical rho for 1/sqrt(2)*(|HHHH>+ i |VVVV>)
% rho_tr([1 16],[1 16])=[0.5 0.5i;-0.5i 0.5]; 
% %2. Theoretical rho for loop cluster
qinit = [1 0]';
Hadamad =   [1 1; 1 -1]/sqrt(2);
q_plus = qo(Hadamad*qinit);
graph_ini = tensor(tensor(q_plus,q_plus),tensor(q_plus,q_plus));
Q1 = qgate(4,cphase,[1,2]);
Q2 = qgate(4,cphase,[2,4]);
Q4 = qgate(4,cphase,[1,3]);
Q3 = qgate(4,cphase,[1,4]);
% Q3 = qgate(4,cphase,[2,3]);

lpGraph = Q4*Q3*Q2*Q1*graph_ini;
rho_graph = double(lpGraph)*double(lpGraph)';
for ind1=1:256
        n_t(ind1)=100*proj4q(ind1,:)*rho_graph*proj4q(ind1,:)';  %assuming there are 100 4-fold coincidences.
end
disp('Loading data...');
folder = [p1(1:strfind(p1, 'Dropbox')-1) 'Dropbox\Data\Tomo\2011-08-31\USE\'];
%folder = 'C:\Dropbox\Data\Tomo\2011-08-31\USE\';
dirListing = dir([folder '*.txt']);
for d = 1:length(dirListing)
    fileName = fullfile(folder,dirListing(d).name); 
    raw_counts = load(fileName);
    coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.92;
    coinc(find(coinc<0))=0.1;
    
    if size(findstr(dirListing(d).name,'HH'),2) == 1 %pol_proj(1), n(pol_proj*16-16+(1:16))
        n(1*16-16+(1:16)') = coinc;
        errload = errload + 1;
    end
    if size(findstr(dirListing(d).name,'HV'),2) == 1  %pol_proj(2)
        n(2*16-16+(1:16)') = coinc;
        errload = errload + 1;        
    end
    if size(findstr(dirListing(d).name,'VV'),2) == 1  %pol_proj(3)
        n(3*16-16+(1:16)') = coinc;
        errload = errload + 1;        
    end
    if size(findstr(dirListing(d).name,'VH'),2) == 1  %pol_proj(4)
        n(4*16-16+(1:16)') = coinc;
        errload = errload + 1;        
    end
    if size(findstr(dirListing(d).name,'RH'),2) == 1  %pol_proj(5)
        n(5*16-16+(1:16)') = coinc;
        errload = errload + 1;        
    end
    if size(findstr(dirListing(d).name,'RV'),2) == 1  %pol_proj(6)
        n(6*16-16+(1:16)') = coinc;
        errload = errload + 1;        
    end
    if size(findstr(dirListing(d).name,'DV'),2) == 1  %pol_proj(7)
        n(7*16-16+(1:16)') = coinc;
        errload = errload + 1;        
    end
    if size(findstr(dirListing(d).name,'DH'),2) == 1  %pol_proj(8)
        n(8*16-16+(1:16)') = coinc;
        errload = errload + 1;        
    end
    if size(findstr(dirListing(d).name,'DR'),2) == 1  %pol_proj(9)
        n(9*16-16+(1:16)') = coinc;
        errload = errload + 1;        
    end
    if size(findstr(dirListing(d).name,'DD'),2) == 1  %pol_proj(10)
        n(10*16-16+(1:16)') = coinc;
        errload = errload + 1;        
    end
    if size(findstr(dirListing(d).name,'RD'),2) == 1  %pol_proj(11)
        n(11*16-16+(1:16)') = coinc;
        errload = errload + 1;        
    end
    if size(findstr(dirListing(d).name,'HD'),2) == 1  %pol_proj(12)
        n(12*16-16+(1:16)') = coinc;
        errload = errload + 1;        
    end
    if size(findstr(dirListing(d).name,'VD'),2) == 1  %pol_proj(13)
        n(13*16-16+(1:16)') = coinc;
        errload = errload + 1;        
    end
    if size(findstr(dirListing(d).name,'VL'),2) == 1  %pol_proj(14)
        n(14*16-16+(1:16)') = coinc;
        errload = errload + 1;        
    end
    if size(findstr(dirListing(d).name,'HL'),2) == 1  %pol_proj(15)
        n(15*16-16+(1:16)') = coinc;
        errload = errload + 1;        
    end
    if size(findstr(dirListing(d).name,'RL'),2) == 1  %pol_proj(16)
        n(16*16-16+(1:16)') = coinc;
        errload = errload + 1;        
    end   
end

if errload < 16
    errload
    error('Data loading not successful. Some basis missing.');        
else if errload >16
        error('Data loading not successful. Multiple files found in the same basis.');
    else if find(n==0)
            error('16 Basis not fully loaded. Some basis missing while others with multiples');
        end
    end
end
ndisp = reshape(n,16,16);

disp('Starting linear inverstion...');
for nu=1:1:256
    rho = rho + reshape(M(nu,:,:),16,16)*n(nu);
    rho_tr = rho_tr + reshape(M(nu,:,:),16,16)*n_t(nu);
end

rho = rho/sum(n([1:4,17:20,33:36,49:52])); %Renomalization, in H, V basis. HHHH, HHHV, etc.
rho_tr = rho_tr/sum(n_t([1:4,17:20,33:36,49:52]));


figure(1);
subplot(2,2,1);
bar3(real(rho));
title('real part - exp, LI');

subplot(2,2,2);
bar3(imag(rho));
title('imaginary part - exp, LI');

subplot(2,2,3);
bar3(real(rho_tr));
title('real part - theory');

subplot(2,2,4);
bar3(imag(rho_tr));
title('imaginary part - theory');


d1=eig(rho);
purity1 = trace(rho*rho);

fidelity1 = fidelity(rho_tr,rho);

disp('Starting MLH using SeDumi, initializing...');
%Maximum likelyhood estimation, using SeDuMi.
oper4q=zeros(256,16,16);
for ind=1:256
    oper4q(ind,:,:) = proj4q(ind,:)'*proj4q(ind,:);
end
if MLHFlag
    rhosdm=sdpvar(16,16,'hermitian','complex');
    n_sdm=sdpvar(256,1);
    assign(rhosdm,eye(16));
    norm1 = sum(n(1:4)+n(16+(1:4))+n(32+(1:4))+n(48+(1:4)));
    val=0;
    for ind=1:256
        n_sdm(ind)=trace(reshape(oper4q(ind,:,:),16,16)*rhosdm);
    end
    val=sum(n_sdm-2*n/norm1+(n/norm1).^2./double(n_sdm));
    %val=sum((n_sdm-n/norm1).^2);

%     
%     t0=FindInitialT4q(rho);
%     tsdm=sdpvar(16,16,'full','complex');
%     n_sdm=sdpvar(256,1);
%     assign(tsdm,t0);
%     norm1 = sum(n(1:4)+n(16+(1:4))+n(32+(1:4))+n(48+(1:4)));
%     val=0;
%     rhosdm=tsdm*tsdm';
%     for ind=1:256
%         n_sdm(ind)=trace(reshape(oper4q(ind,:,:),16,16)*rhosdm);
%     end
%     val=sum(n_sdm-2*n/norm1+(n/norm1).^2./double(n_sdm));
    
    disp('Start the solver...');
    solvesdp([rhosdm>0; trace(rhosdm)==1],val);
    rho_sdm = double(rhosdm);
    rho_sdm=rho_sdm/trace(rho_sdm);

    d2=eig(rho_sdm);
    purity2 = trace(rho_sdm*rho_sdm);
    
    figure(2); 
    subplot(1,2,1);bar3(real(rho_sdm));
    title('Real (MLH Estimation_{SDM})');
    subplot(1,2,2);bar3(imag(rho_sdm));
    title('Imag (MLH Estimation_{SDM})');
    disp(['Eig_LI','   Eig_MLH']);
    disp(real([d1,d2]));
    disp(['Purity_LI','  Purity_MLH']);
    disp(real([purity1,purity2]));
    disp(['Fidelity_LI_Th ' 'SDM_Th ' 'LI_SDM']);
    disp(real([fidelity(rho,rho_tr) fidelity(rho_sdm,rho_tr)...
        fidelity(rho_sdm,rho)]));
else
    disp('Eig_LI');
    disp(real(d1));
    disp('Purity_LI');
    disp(real(purity1));
    disp('Fidelity_LI');
    disp(real(fidelity1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



norm_nt=sum(n_t([1:4,17:20,33:36,49:52]));
norm_n =sum(n([1:4,17:20,33:36,49:52]));
err1 = (n_t'/norm_nt-n/norm_n)/16; %normalized error in probability.
figure(4);plot(err1);
ind=find(abs(err1)>0.007);
r_ind=rem(ind,16);
figure(5);hist(r_ind,max(r_ind)-min(r_ind));
