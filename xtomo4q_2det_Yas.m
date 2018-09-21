MLHFlag=0; %Maximum likelihood not implemented yet.

proj4q=zeros(256,16);
M = zeros(256,16,16);
global n;
n = zeros(256,1);
rho = zeros(16,16);
counting_time = 15;
coin_win = 190e-9;
errload = 0;

%%%% Construct the 4-qubit projectors from 2-qubit projectors. This is a
%%%% bit time consuming, but it only need to be done once.
proj1=proj_2q_pol;
proj2 = proj_2q_path('uu');
proj2=proj_2q_pol;
for ind1=1:16
    for ind2=1:16
            proj4q(ind1*16-16+ind2,:)=tensor_product(proj1(ind1,:),proj2(ind2,:));
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

% qinit = [1 0]';
% Hadamad =   [1 1; 1 -1]/sqrt(2);
% %q_plus = qo(Hadamad*qinit);
% q_plus = Hadamad*qinit;
% graph_ini = kron(kron(q_plus,q_plus),kron(q_plus,q_plus));
% Q1 = qgate(4,cphase,[1,2]);
% Q2 = qgate(4,cphase,[2,4]);
% Q3 = qgate(4,cphase,[1,4]);
% Q4 = qgate(4,cphase,[1,3]);
% lpGraph = Q4*Q3*Q2*Q1*graph_ini;
% lpGraph = 1/4*[1, 1, 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, -1, -1, 1, 1];
% rho_t = double(lpGraph)*double(lpGraph)';
% for ind1=1:256
%        n_t(ind1)=100*proj4q(ind1,:)*rho_t*proj4q(ind1,:)';  %assuming there are 100 4-fold coincidences.
% end
% n=n_t;

%p1 = pwd;
% folder = load('C:\Documents and Settings\Jazz\My Documents\Thesis\Use for Tomo\');
% dirListing = dir([folder '2011-08-07_afterInterfero*']);
% for d = 1:length(dirListing)
%     fileName = fullfile(folder,dirListing(d).name); 
%     raw_counts = load(fileName);
%     coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.96;
%     coinc(find(coinc<0))=0.01;
%     
%     if size(findstr(dirListing(d).name,'HH'),2) == 1 %pol_proj(1), n(pol_proj*16-16+(1:16))
%         n(1*16-16+(1:16)') = coinc;
%         errload = errload + 1;
%     end
%     if size(findstr(dirListing(d).name,'HV'),2) == 1  %pol_proj(2)
%         n(2*16-16+(1:16)') = coinc;
%         errload = errload + 1;        
%     end
%     if size(findstr(dirListing(d).name,'VV'),2) == 1  %pol_proj(3)
%         n(3*16-16+(1:16)') = coinc;
%         errload = errload + 1;        
%     end
%     if size(findstr(dirListing(d).name,'VH'),2) == 1  %pol_proj(4)
%         n(4*16-16+(1:16)') = coinc;
%         errload = errload + 1;        
%     end
%     if size(findstr(dirListing(d).name,'RH'),2) == 1  %pol_proj(5)
%         n(5*16-16+(1:16)') = coinc;
%         errload = errload + 1;        
%     end
%     if size(findstr(dirListing(d).name,'RV'),2) == 1  %pol_proj(6)
%         n(6*16-16+(1:16)') = coinc;
%         errload = errload + 1;        
%     end
%     if size(findstr(dirListing(d).name,'DV'),2) == 1  %pol_proj(7)
%         n(7*16-16+(1:16)') = coinc;
%         errload = errload + 1;        
%     end
%     if size(findstr(dirListing(d).name,'DH'),2) == 1  %pol_proj(8)
%         n(8*16-16+(1:16)') = coinc;
%         errload = errload + 1;        
%     end
%     if size(findstr(dirListing(d).name,'DR'),2) == 1  %pol_proj(9)
%         n(9*16-16+(1:16)') = coinc;
%         errload = errload + 1;        
%     end
%     if size(findstr(dirListing(d).name,'RR'),2) == 1  %pol_proj(10)
%         n(10*16-16+(1:16)') = coinc;
%         errload = errload + 1;        
%     end
%     if size(findstr(dirListing(d).name,'RD'),2) == 1  %pol_proj(11)
%         n(11*16-16+(1:16)') = coinc;
%         errload = errload + 1;        
%     end
%     if size(findstr(dirListing(d).name,'HD'),2) == 1  %pol_proj(12)
%         n(12*16-16+(1:16)') = coinc;
%         errload = errload + 1;        
%     end
%     if size(findstr(dirListing(d).name,'VD'),2) == 1  %pol_proj(13)
%         n(13*16-16+(1:16)') = coinc;
%         errload = errload + 1;        
%     end
%     if size(findstr(dirListing(d).name,'VL'),2) == 1  %pol_proj(14)
%         n(14*16-16+(1:16)') = coinc;
%         errload = errload + 1;        
%     end
%     if size(findstr(dirListing(d).name,'HL'),2) == 1  %pol_proj(15)
%         n(15*16-16+(1:16)') = coinc;
%         errload = errload + 1;        
%     end
%     if size(findstr(dirListing(d).name,'RL'),2) == 1  %pol_proj(16)
%         n(16*16-16+(1:16)') = coinc;
%         errload = errload + 1;        
%     end   
% end
% 
% if errload < 16
%     error('Data loading not successful. Some basis missing.');
% else if errload >16
%         error('Data loading not successful. Multiple files found in the same basis.');
%     else if find(n==0)
%             error('16 Basis not fully loaded. Some basis missing while others with multiples');
%         end
%     end
% end

raw_counts =load('C:\Documents and Settings\Jazz\My Documents\Thesis\Use for Tomo\2011-08-07_afterInterfero_corrected_maxDD_HH_15s_01.txt')
coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.96;
    coinc(find(coinc<0))=0.01;
    n(1*16-16+(1:16)') = coinc;
raw_counts =load('C:\Documents and Settings\Jazz\My Documents\Thesis\Use for Tomo\2011-08-07_afterInterfero_corrected_maxDD_HV_15s_01.txt')
coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.96;
    coinc(find(coinc<0))=0.01;
    n(2*16-16+(1:16)') = coinc;

raw_counts =load('C:\Documents and Settings\Jazz\My Documents\Thesis\Use for Tomo\2011-08-07_afterInterfero_corrected_maxDD_VH_15s_01.txt')
coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.96;
    coinc(find(coinc<0))=0.01;
    n(3*16-16+(1:16)') = coinc;
raw_counts =load('C:\Documents and Settings\Jazz\My Documents\Thesis\Use for Tomo\2011-08-07_afterInterfero_corrected_maxDD_VV_15s_01.txt')
coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.96;
    coinc(find(coinc<0))=0.01;
    n(4*16-16+(1:16)') = coinc;
raw_counts =load('C:\Documents and Settings\Jazz\My Documents\Thesis\Use for Tomo\2011-08-07_afterInterfero_corrected_maxDD_RH_15s_01.txt')
coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.96;
    coinc(find(coinc<0))=0.01;
    n(5*16-16+(1:16)') = coinc;
raw_counts =load('C:\Documents and Settings\Jazz\My Documents\Thesis\Use for Tomo\2011-08-07_afterInterfero_corrected_maxDD_RV_15s_01.txt')
coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.96;
    coinc(find(coinc<0))=0.01;
    n(6*16-16+(1:16)') = coinc;
raw_counts =load('C:\Documents and Settings\Jazz\My Documents\Thesis\Use for Tomo\2011-08-07_afterInterfero_corrected_maxDD_DV_15s_01.txt')
coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.96;
    coinc(find(coinc<0))=0.01;
    n(7*16-16+(1:16)') = coinc;
raw_counts =load('C:\Documents and Settings\Jazz\My Documents\Thesis\Use for Tomo\2011-08-07_afterInterfero_corrected_maxDD_DH_15s_01.txt')
coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.96;
    coinc(find(coinc<0))=0.01;
    n(8*16-16+(1:16)') = coinc;
raw_counts =load('C:\Documents and Settings\Jazz\My Documents\Thesis\Use for Tomo\2011-08-07_afterInterfero_corrected_maxDD_DR_15s_01.txt')
coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.96;
    coinc(find(coinc<0))=0.01;
    n(9*16-16+(1:16)') = coinc;
raw_counts =load('C:\Documents and Settings\Jazz\My Documents\Thesis\Use for Tomo\2011-08-07_afterInterfero_corrected_maxDD_RR_15s_01.txt')
coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.96;
    coinc(find(coinc<0))=0.01;
    n(10*16-16+(1:16)') = coinc;
raw_counts =load('C:\Documents and Settings\Jazz\My Documents\Thesis\Use for Tomo\2011-08-07_afterInterfero_corrected_maxDD_RD_15s_01.txt')
coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.96;
    coinc(find(coinc<0))=0.01;
    n(11*16-16+(1:16)') = coinc;
raw_counts =load('C:\Documents and Settings\Jazz\My Documents\Thesis\Use for Tomo\2011-08-07_afterInterfero_corrected_maxDD_HD_15s_01.txt')
coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.96;
    coinc(find(coinc<0))=0.01;
    n(12*16-16+(1:16)') = coinc;
raw_counts =load('C:\Documents and Settings\Jazz\My Documents\Thesis\Use for Tomo\2011-08-07_afterInterfero_corrected_maxDD_VD_15s_01.txt')
coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.96;
    coinc(find(coinc<0))=0.01;
    n(13*16-16+(1:16)') = coinc;
raw_counts =load('C:\Documents and Settings\Jazz\My Documents\Thesis\Use for Tomo\2011-08-07_afterInterfero_corrected_maxDD_VL_15s_01.txt')
coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.96;
    coinc(find(coinc<0))=0.01;
    n(14*16-16+(1:16)') = coinc;
raw_counts =load('C:\Documents and Settings\Jazz\My Documents\Thesis\Use for Tomo\2011-08-07_afterInterfero_corrected_maxDD_HL_15s_01.txt')
coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.96;
    coinc(find(coinc<0))=0.01;
    n(15*16-16+(1:16)') = coinc;
raw_counts =load('C:\Documents and Settings\Jazz\My Documents\Thesis\Use for Tomo\2011-08-07_afterInterfero_corrected_maxDD_RL_15s_01.txt')
coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.96;
    coinc(find(coinc<0))=0.01;
    n(16*16-16+(1:16)') = coinc;

ndisp = reshape(n,16,16)

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




